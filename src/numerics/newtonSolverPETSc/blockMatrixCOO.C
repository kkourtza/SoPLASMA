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
    nFields_(nFields)
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
    const label rowOff = rowField*nCells_;
    const label colOff = colField*nCells_;

    forAll(coeff, c)
    {
        add(rowOff + c, colOff + c, coeff[c]*scaling);
    }
}


void Foam::blockMatrixCOO::addFvMatrix
(
    const label field,
    const fvScalarMatrix& m,
    const scalar scaling
)
{
    if (Pstream::parRun())
    {
        FatalErrorInFunction
            << "Assembling a block Jacobian in parallel is not implemented."
            << nl << nl
            << "    Processor-interface faces contribute OFF-DIAGONAL entries"
            << " that couple cells owned by different ranks, and they are not"
            << " gathered here. Assembling without them would silently"
            << " produce a matrix describing a set of DISCONNECTED"
            << " subdomains -- a preconditioner built on it would look"
            << " plausible and be wrong, which is worse than refusing." << nl
            << "    petsc4Foam's own buildMat() handles these via"
            << " lduAddressing::patchAddr() and the interface"
            << " boundaryCoeffs; follow that when adding parallel support."
            << nl << exit(FatalError);
    }

    const lduAddressing& addr = m.lduAddr();
    const labelUList& upp = addr.upperAddr();
    const labelUList& low = addr.lowerAddr();

    const label off = field*nCells_;

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
        add(off + c, off + c, d[c]*scaling);
    }

    // OFF-DIAGONAL, both triangles. An asymmetric matrix carries a separate
    // lower(); a symmetric one reuses upper(), which is the same convention
    // petscSolver::buildMat uses.
    const scalarField& uppVal = m.upper();
    const scalarField& lowVal = (m.hasLower() ? m.lower() : m.upper());

    forAll(upp, f)
    {
        // row = owner (lower-numbered cell), col = neighbour
        add(off + low[f], off + upp[f], uppVal[f]*scaling);

        // and the transpose position
        add(off + upp[f], off + low[f], lowVal[f]*scaling);
    }
}


// ************************************************************************* //
