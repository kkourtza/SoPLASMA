/*---------------------------------------------------------------------------*\
License
    This file is part of SoPLASMA.

    Copyright (C) 2026

    This program is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

Class
    Foam::blockMatrixCOO

Description
    Accumulates a multi-field block Jacobian in COO (coordinate) form, ready to
    be handed across the plain-types boundary in `snesBridge.H` and turned into
    a PETSc MATAIJ there.

    WHY THIS EXISTS AT ALL. The physics-based preconditioner assembles each
    equation's linearisation with `fvm::` operators already; what was missing
    was any representation of the CROSS-COUPLING between fields, and PETSc's
    PCFIELDSPLIT -- the thing that can exploit it -- cannot work from a
    matrix-free operator. It needs a genuinely assembled preconditioning
    matrix. See docs/design/newton-outer-solver-design.md.

    WHY NOT petsc4Foam. `petscSolver::buildMat()` is exactly the ldu2csr
    conversion we want and it is exported, but it sizes the matrix to ONE
    field (`MatSetSizes(A, nrows, nrows, ...)` with `nrows = lduAddr.size()`),
    with no field-block offset, so it cannot express a multi-field system. Its
    COO LAYOUT is followed here exactly, including the row/column offset it
    already uses for multi-region coupled assembly -- fields here in place of
    regions.

    THE CONVENTION THAT IS EASY TO GET WRONG. petsc4Foam copies
    `matrix.diag()` directly, which is right only because
    `fvMatrix::solveSegregated()` has already folded the boundary contribution
    into the diagonal via `addBoundaryDiag()` before the lduMatrix solver ever
    sees it. Converting an `fvScalarMatrix` ourselves, we must do that
    explicitly on a COPY of the diagonal, or every boundary cell is wrong.

    SCALING. Everything PETSc sees in this solver is scaled (`x_p = x/sX`,
    `F_p = F/sF`), so the Jacobian must be scaled to match: an entry in block
    (row field r, column field c) is multiplied by `sX[c]/sF[r]`. A Pmat
    without this is a perfectly valid matrix of the WRONG system.

SourceFiles
    blockMatrixCOO.C

\*---------------------------------------------------------------------------*/

#ifndef blockMatrixCOO_H
#define blockMatrixCOO_H

#include "fvMatrices.H"
#include "DynamicList.H"
#include "globalIndex.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class blockMatrixCOO Declaration
\*---------------------------------------------------------------------------*/

class blockMatrixCOO
{
    //- Local cell count -- the size of ONE field block.
    const label nCells_;

    //- Number of field blocks.
    const label nFields_;

    //- RAGGED TAIL: extra rows appended AFTER the nFields uniform blocks,
    //  for cells that exist for ONE field only.
    //
    //  multiRegionPoisson solves ePotential over the gas AND every dielectric
    //  or farField region, while the species and energy exist in the GAS
    //  alone. Padding every block to the total was rejected: a `farField`
    //  region is "a FICTITIOUS air/vacuum region: Poisson only" and is
    //  EXPECTED to be larger than the gas ("Poisson is long-range and usually
    //  needs a larger domain than the plasma does"), so padding would multiply
    //  the DOF count several-fold for nothing.
    //
    //  A tail costs exactly the non-gas cells and leaves every existing index
    //  expression -- field*nCells + cell -- untouched. It is ZERO for a
    //  single-region case, which is what makes those bit-identical.
    //
    //      [ f0 gas ][ f1 gas ] ... [ fN gas ][ f0 non-gas ]
    //                                          ^ nTail_ rows
    const label nTail_;

    //- GLOBAL row numbering across ranks, over this rank's nFields*nCells
    //  rows. PETSc's Vec/Mat layout is rank-contiguous, and our per-rank
    //  layout is field-major inside that block, so the global row of
    //  (field f, local cell c) is localStart() + f*nCells + c. In serial
    //  localStart() is 0 and every index below is unchanged.
    const globalIndex rowNumbering_;

    DynamicList<label> rows_;
    DynamicList<label> cols_;
    DynamicList<scalar> vals_;


public:

    //- Construct for a given block layout. nTail defaults to 0, so every
    //  existing caller and every single-region case is unchanged.
    blockMatrixCOO
    (
        const label nCells,
        const label nFields,
        const label nTail = 0
    );

    //- Discard all accumulated entries, keeping the allocation
    void clear();

    //- Total unknowns, i.e. the matrix dimension
    label nRows() const noexcept { return nFields_*nCells_ + nTail_; }

    //- Size of the ragged tail (0 when single-region)
    label nTail() const noexcept { return nTail_; }

    label nEntries() const noexcept { return rows_.size(); }

    const label* rows() const noexcept { return rows_.cdata(); }
    const label* cols() const noexcept { return cols_.cdata(); }
    const scalar* vals() const noexcept { return vals_.cdata(); }

    //- Add one entry, in GLOBAL (block-offset) indices. Duplicate
    //  coordinates are summed by PETSc's COO assembly, which is what makes
    //  contributions from separate blocks composable.
    void add(const label row, const label col, const scalar v);

    //- GLOBAL row/column index of (field, local cell)
    label globalRow(const label field, const label cell) const
    {
        return rowNumbering_.localStart() + field*nCells_ + cell;
    }

    //- GLOBAL row/column index of a TAIL cell (a non-gas cell of field 0).
    //  `tailCell` is 0-based within the tail, i.e. assembledCell - nCells_
    //  using the assembly's own cellOffsets().
    label globalTailRow(const label tailCell) const
    {
        return rowNumbering_.localStart() + nFields_*nCells_ + tailCell;
    }

    //- GLOBAL row of an ASSEMBLED cell index for field 0 (ePotential):
    //  gas cells map into block 0, non-gas cells into the tail. This is the
    //  one mapping the multi-region Poisson block needs, and it is the
    //  IDENTITY for a single-region case.
    label globalPoissonRow(const label assembledCell) const
    {
        return (assembledCell < nCells_)
             ? globalRow(0, assembledCell)
             : globalTailRow(assembledCell - nCells_);
    }

    //- Add a DIAGONAL coupling block: d(residual of field rowField)/d(field
    //  colField), given cell-wise. Used for the Poisson<->species coupling,
    //  which is diagonal because chargeDensity is a pointwise sum of the
    //  species densities.
    void addDiagonalBlock
    (
        const label rowField,
        const label colField,
        const scalarField& coeff,
        const scalar scaling
    );

    //- Add a full fvScalarMatrix as the (field, field) diagonal block.
    //
    //  Follows petsc4Foam's COO layout: the diagonal, then one entry per
    //  internal face for the upper triangle and one for the lower. The
    //  diagonal has addBoundaryDiag() applied first, exactly as
    //  fvMatrix::solveSegregated() does before handing a matrix to a solver.
    //
    //  Processor interfaces ARE contributed: each coupled patch adds one
    //  off-diagonal entry per face, at (this cell, the NEIGHBOUR's global
    //  row), with value -boundaryCoeffs. The neighbour's global row arrives
    //  through OpenFOAM's own interface transfer rather than any manual rank
    //  arithmetic, exactly as petsc4Foam's buildMat does it.
    //  rowFactor multiplies every entry of row c. Pass 1/V to convert an
    //  fvm:: matrix (which is VOLUME-INTEGRATED -- e.g. EulerDdtScheme sets
    //  diag = rDeltaT*mesh().V()) into the per-unit-volume convention the
    //  fvc::-built residual uses. Mixing the two silently yields a
    //  preconditioner wrong by a per-row factor of V, which is invisible on a
    //  uniform mesh and wrong on a graded one.
    void addFvMatrix
    (
        const label field,
        const fvScalarMatrix& m,
        const scalar scaling,
        const scalarField& rowFactor
    );

    //- Place a Poisson-only region's matrix into the RAGGED TAIL rows.
    //
    //  tailOffset is the region's first row COUNTED FROM THE START OF THE
    //  TAIL, so consecutive regions pass their running cell offset.
    //
    //  The block is DIAGONAL in the tail: the inter-region coupling is carried
    //  by the explicit residual (coupledElectricPotential's Robin condition),
    //  not by this matrix, which exists only to give the preconditioner a real
    //  elliptic operator on those rows. WITHOUT IT THOSE ROWS ARE EMPTY, the
    //  phi split's diagonal is zero there, and the first PCApply raises an FPE
    //  -- measured on needleDBD 2026-09-11, after the residual itself had
    //  already evaluated cleanly (SNES Function norm 1.037e+03).
    void addFvMatrixTail
    (
        const label tailOffset,
        const fvScalarMatrix& m,
        const scalar scaling,
        const scalarField& rowFactor
    );

    //- Place an fvMatrix in an OFF-DIAGONAL block: rows belong to rowField,
    //  columns to colField.
    //
    //  This is what lets the Pmat carry d(species residual)/d(ePotential) --
    //  the response of the drift flux to the field. Without it the Pmat is
    //  block-triangular in {phi, transport}, and PCFIELDSPLIT's Schur
    //  complement S = A_tt - A_tp*A_pp^-1*A_pt degenerates to A_tt exactly,
    //  because A_tp is zero. The Schur complement then contributes NOTHING,
    //  which defeats the entire reason an assembled Pmat exists here.
    void addFvMatrixBlock
    (
        const label rowField,
        const label colField,
        const fvScalarMatrix& m,
        const scalar scaling,
        const scalarField& rowFactor
    );

    //- PER-CELL scaled variant: the row factor and the COLUMN factor are both
    //  arrays, so a state scaling that varies cell by cell can be represented.
    //
    //  WHY A SEPARATE OVERLOAD rather than generalising the one above. Under a
    //  per-cell state scaling xhat_j = x_j/sX_j and residual scaling
    //  Fhat_i = F_i/sF_i, the Jacobian transforms as
    //      Jhat_ij = J_ij * sX_j / sF_i
    //  so the entry needs the COLUMN's sX and the ROW's sF -- and the column
    //  index differs from the row index on every off-diagonal. The scalar
    //  `scaling` argument of the overload above cannot express that.
    //
    //  It is additive because folding sF into rowFactor and sX into a new
    //  colFactor would CHANGE THE MULTIPLICATION ORDER of the existing path,
    //  and floating-point multiplication is not associative -- the uniform
    //  case would stop being bit-identical to what it is today, which is the
    //  regression gate the per-cell switch is verified against. Leaving the
    //  original untouched makes "switch off changes nothing" provable rather
    //  than argued.
    //
    //  rowFactor[i] must already carry 1/V_i and 1/sF_i; colFactor carries
    //  sX_j. No scalar multiplier: everything is in the two.
    //
    //  colFactor IS A volScalarField, not a scalarField, and that is
    //  load-bearing. A processor-interface entry sits at (local row, NEIGHBOUR
    //  column), so it needs the neighbour's sX -- which lives on another rank.
    //  A volScalarField already carries its processor-patch values, exchanged
    //  by correctBoundaryConditions(), so `colFactor.boundaryField()[patchi]`
    //  IS the neighbour's scaling with no hand-rolled transfer. (lduInterface's
    //  internalFieldTransfer only handles labels, so the alternative would have
    //  been a manual processorLduInterface send/receive -- working inside
    //  OpenFOAM's own structures instead is rule 40.)
    //
    //  The CALLER must have called correctBoundaryConditions() on it.
    void addFvMatrixBlockPerCell
    (
        const label rowField,
        const label colField,
        const fvScalarMatrix& m,
        const scalarField& rowFactor,
        const volScalarField& colFactor
    );
};

} // End namespace Foam

#endif

// ************************************************************************* //
