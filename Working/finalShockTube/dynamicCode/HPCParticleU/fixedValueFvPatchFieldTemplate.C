/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) YEAR OpenFOAM Foundation
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

#include "fixedValueFvPatchFieldTemplate.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "unitConversion.H"
//{{{ begin codeInclude

//}}} end codeInclude


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

extern "C"
{
    // dynamicCode:
    // SHA1 = cb97bb81a9b2557ddbcbad5e4e7952fe79b923b6
    //
    // unique function name that can be checked if the correct library version
    // has been loaded
    void HPCParticleU_cb97bb81a9b2557ddbcbad5e4e7952fe79b923b6(bool load)
    {
        if (load)
        {
            // code that can be explicitly executed after loading
        }
        else
        {
            // code that can be explicitly executed before unloading
        }
    }
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

makeRemovablePatchTypeField
(
    fvPatchVectorField,
    HPCParticleUFixedValueFvPatchVectorField
);


const char* const HPCParticleUFixedValueFvPatchVectorField::SHA1sum =
    "cb97bb81a9b2557ddbcbad5e4e7952fe79b923b6";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

HPCParticleUFixedValueFvPatchVectorField::
HPCParticleUFixedValueFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF)
{
    if (false)
    {
        Info<<"construct HPCParticleU sha1: cb97bb81a9b2557ddbcbad5e4e7952fe79b923b6"
            " from patch/DimensionedField\n";
    }
}


HPCParticleUFixedValueFvPatchVectorField::
HPCParticleUFixedValueFvPatchVectorField
(
    const HPCParticleUFixedValueFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper)
{
    if (false)
    {
        Info<<"construct HPCParticleU sha1: cb97bb81a9b2557ddbcbad5e4e7952fe79b923b6"
            " from patch/DimensionedField/mapper\n";
    }
}


HPCParticleUFixedValueFvPatchVectorField::
HPCParticleUFixedValueFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict)
{
    if (false)
    {
        Info<<"construct HPCParticleU sha1: cb97bb81a9b2557ddbcbad5e4e7952fe79b923b6"
            " from patch/dictionary\n";
    }
}


HPCParticleUFixedValueFvPatchVectorField::
HPCParticleUFixedValueFvPatchVectorField
(
    const HPCParticleUFixedValueFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf)
{
    if (false)
    {
        Info<<"construct HPCParticleU sha1: cb97bb81a9b2557ddbcbad5e4e7952fe79b923b6"
            " as copy\n";
    }
}


HPCParticleUFixedValueFvPatchVectorField::
HPCParticleUFixedValueFvPatchVectorField
(
    const HPCParticleUFixedValueFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF)
{
    if (false)
    {
        Info<<"construct HPCParticleU sha1: cb97bb81a9b2557ddbcbad5e4e7952fe79b923b6 "
            "as copy/DimensionedField\n";
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

HPCParticleUFixedValueFvPatchVectorField::
~HPCParticleUFixedValueFvPatchVectorField()
{
    if (false)
    {
        Info<<"destroy HPCParticleU sha1: cb97bb81a9b2557ddbcbad5e4e7952fe79b923b6\n";
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void HPCParticleUFixedValueFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (false)
    {
        Info<<"updateCoeffs HPCParticleU sha1: cb97bb81a9b2557ddbcbad5e4e7952fe79b923b6\n";
    }

//{{{ begin code
    #line 39 "/home/anthonygay1812/OpenFOAM/Working/finalShockTube/system/particle/codeDict.HPCParticleU"
const fvMesh& fluidMeshRef = db().parent().lookupObject<fvMesh>("fluid");
    const scalar idx=fluidMeshRef.nCells()-1;

    const volVectorField& UFluid = fluidMeshRef.lookupObject<volVectorField>("U");     
    Info<< nl << "DSMC U BC";
    Info<< nl << "Calculated U: "<< (UFluid[idx]);
    operator==(UFluid[idx]);
//}}} end code

    this->fixedValueFvPatchField<vector>::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

