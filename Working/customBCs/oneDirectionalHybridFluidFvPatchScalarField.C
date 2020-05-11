/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "oneDirectionalHybridFluidFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::oneDirectionalHybridFluidFvPatchScalarField::oneDirectionalHybridFluidFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    fluidIdx_(0),
    particleIdx_(1),
    w_(0)
{}


Foam::oneDirectionalHybridFluidFvPatchScalarField::oneDirectionalHybridFluidFvPatchScalarField
(
    const oneDirectionalHybridFluidFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    fluidIdx_(ptf.fluidIdx_),
    particleIdx_(ptf.particleIdx_),
    w_(ptf.w_)
{}


Foam::oneDirectionalHybridFluidFvPatchScalarField::oneDirectionalHybridFluidFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict, false),
    fluidIdx_(readLabel(dict.lookup("fluidIndex"))),
    particleIdx_(readLabel(dict.lookup("particleIndex"))),
    w_(readScalar(dict.lookup("w")))

{
    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
}


Foam::oneDirectionalHybridFluidFvPatchScalarField::oneDirectionalHybridFluidFvPatchScalarField
(
    const oneDirectionalHybridFluidFvPatchScalarField& tppsf
)
:
    fixedValueFvPatchScalarField(tppsf),
    fluidIdx_(tppsf.fluidIdx_),
    particleIdx_(tppsf.particleIdx_),
    w_(tppsf.w_)
{}


Foam::oneDirectionalHybridFluidFvPatchScalarField::oneDirectionalHybridFluidFvPatchScalarField
(
    const oneDirectionalHybridFluidFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF),
    fluidIdx_(tppsf.fluidIdx_),
    particleIdx_(tppsf.particleIdx_),
    w_(tppsf.w_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::oneDirectionalHybridFluidFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchScalarField::autoMap(m);
}


void Foam::oneDirectionalHybridFluidFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);
}


void Foam::oneDirectionalHybridFluidFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    //const fvPatchVectorField& Up =
    //    patch().lookupPatchField<volVectorField, vector>(UName_);
    //const fvMesh& fluidMeshRef_ = patch().boundaryMesh().mesh();
    const volScalarField& pFluid_ = db().lookupObject<volScalarField>("p"); 
    const volScalarField& pParticle_ = db().lookupObject<volScalarField>("p_Mean");
    
    const volScalarField& pFluidOld_=pFluid_.oldTime();

    operator==
    (
        ((1-w_)*pFluidOld_[fluidIdx_]+w_*(0.5*pFluid_[fluidIdx_]+0.5*pParticle_[particleIdx_]))
    );

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::oneDirectionalHybridFluidFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
   // writeEntry(os, "w", w_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        oneDirectionalHybridFluidFvPatchScalarField
    );
}

// ************************************************************************* //
