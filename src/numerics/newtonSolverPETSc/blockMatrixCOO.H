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
    matrix. See doc/newton-outer-solver-design.md.

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

    //- Construct for a given block layout
    blockMatrixCOO(const label nCells, const label nFields);

    //- Discard all accumulated entries, keeping the allocation
    void clear();

    //- Total unknowns, i.e. the matrix dimension
    label nRows() const noexcept { return nFields_*nCells_; }

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
};

} // End namespace Foam

#endif

// ************************************************************************* //
