/*---------------------------------------------------------------------------*\
License
    This file is part of SoPLASMA.

    Copyright (C) 2026

    This program is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

\*---------------------------------------------------------------------------*/

#include "blockMatrixCOO.H"
#include "lduAddressing.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blockMatrixCOO::blockMatrixCOO
(
    const label nCells,
    const label nFields
)
:
    nCells_(nCells),
    nFields_(nFields),
    rowNumbering_(nFields*nCells)
{
    // Diagonal for every unknown, plus both triangles of every internal face
    // for every field block. An estimate only -- DynamicList grows if the
    // coupling blocks add more.
    const label estimate = nFields*nCells*8;
    rows_.reserve(estimate);
    cols_.reserve(estimate);
    vals_.reserve(estimate);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::blockMatrixCOO::clear()
{
    rows_.clear();
    cols_.clear();
    vals_.clear();
}


void Foam::blockMatrixCOO::add
(
    const label row,
    const label col,
    const scalar v
)
{
    rows_.append(row);
    cols_.append(col);
    vals_.append(v);
}


void Foam::blockMatrixCOO::addDiagonalBlock
(
    const label rowField,
    const label colField,
    const scalarField& coeff,
    const scalar scaling
)
{
    forAll(coeff, c)
    {
        add(globalRow(rowField, c), globalRow(colField, c), coeff[c]*scaling);
    }
}


void Foam::blockMatrixCOO::addFvMatrix
(
    const label field,
    const fvScalarMatrix& m,
    const scalar scaling,
    const scalarField& rowFactor
)
{
    // A diagonal block is the off-diagonal case with rowField == colField.
    addFvMatrixBlock(field, field, m, scaling, rowFactor);
}


void Foam::blockMatrixCOO::addFvMatrixBlock
(
    const label rowField,
    const label colField,
    const fvScalarMatrix& m,
    const scalar scaling,
    const scalarField& rowFactor
)
{
    const lduAddressing& addr = m.lduAddr();
    const labelUList& upp = addr.upperAddr();
    const labelUList& low = addr.lowerAddr();

    // DIAGONAL, with the boundary contribution folded in.
    //
    // This is the step petsc4Foam does not appear to do only because
    // fvMatrix::solveSegregated() has already done it by the time a matrix
    // reaches an lduMatrix solver. Working from an fvScalarMatrix directly, we
    // must do it ourselves -- on a COPY, so the caller's matrix is untouched.
    scalarField d(m.diag());
    {
        // fvMatrix::addBoundaryDiag() is PROTECTED, so its (short) body is
        // reproduced here from public API rather than reached by a cast:
        // each patch's internalCoeffs are added into the diagonal at the
        // patch's own cells. Keeping it explicit also documents the sign.
        const lduAddressing& a = m.lduAddr();
        const FieldField<Field, scalar>& intCoeffs = m.internalCoeffs();

        forAll(intCoeffs, patchi)
        {
            const labelUList& pa = a.patchAddr(patchi);
            const scalarField& pc = intCoeffs[patchi];

            forAll(pa, i)
            {
                d[pa[i]] += pc[i];
            }
        }
    }

    forAll(d, c)
    {
        add(globalRow(rowField, c), globalRow(colField, c), d[c]*scaling*rowFactor[c]);
    }

    // OFF-DIAGONAL, both triangles. An asymmetric matrix carries a separate
    // lower(); a symmetric one reuses upper(), which is the same convention
    // petscSolver::buildMat uses.
    const scalarField& uppVal = m.upper();
    const scalarField& lowVal = (m.hasLower() ? m.lower() : m.upper());

    forAll(upp, f)
    {
        // row = owner (lower-numbered cell), col = neighbour
        add
        (
            globalRow(rowField, low[f]), globalRow(colField, upp[f]),
            uppVal[f]*scaling*rowFactor[low[f]]
        );

        // and the transpose position
        add
        (
            globalRow(rowField, upp[f]), globalRow(colField, low[f]),
            lowVal[f]*scaling*rowFactor[upp[f]]
        );
    }

    // ---- PROCESSOR INTERFACES: the entries that couple cells owned by
    // DIFFERENT ranks. Without them the matrix describes a set of
    // disconnected subdomains -- plausible-looking and wrong.
    //
    // The neighbour's GLOBAL row is obtained through OpenFOAM's own interface
    // transfer rather than by computing rank offsets by hand: we send this
    // rank's global row index for every cell adjacent to the interface, and
    // receive the neighbour's. That is exactly how petsc4Foam's buildMat()
    // does it, and it means the field-major-inside-rank layout needs no
    // special treatment -- the number that arrives is already the right one.
    //
    // The transfer carries the global rows OF THIS FIELD, so it is done once
    // per field rather than once per matrix.
    if (Pstream::parRun())
    {
        const lduInterfacePtrsList interfaces(m.psi().mesh().interfaces());

        labelList globalRowsThisField(nCells_);
        forAll(globalRowsThisField, c)
        {
            globalRowsThisField[c] = globalRow(colField, c);
        }

        const label startOfRequests = UPstream::nRequests();

        forAll(interfaces, patchi)
        {
            if (interfaces.set(patchi))
            {
                interfaces[patchi].initInternalFieldTransfer
                (
                    Pstream::commsTypes::nonBlocking,
                    globalRowsThisField
                );
            }
        }

        UPstream::waitRequests(startOfRequests);

        const FieldField<Field, scalar>& bouCoeffs = m.boundaryCoeffs();

        forAll(interfaces, patchi)
        {
            if (!interfaces.set(patchi)) continue;

            const labelUList& faceCells = addr.patchAddr(patchi);

            const labelField nbrRows
            (
                interfaces[patchi].internalFieldTransfer
                (
                    Pstream::commsTypes::nonBlocking,
                    globalRowsThisField
                )
            );

            const scalarField& bc = bouCoeffs[patchi];

            forAll(faceCells, i)
            {
                // MINUS the boundary coefficient: these are THIS side's
                // coefficients, from discretising our own face rather than
                // the neighbour's reversed one (petsc4Foam states the same
                // convention explicitly at its own interface loop).
                add
                (
                    globalRow(rowField, faceCells[i]),
                    nbrRows[i],
                    -bc[i]*scaling*rowFactor[faceCells[i]]
                );
            }
        }
    }
}


// ************************************************************************* //
