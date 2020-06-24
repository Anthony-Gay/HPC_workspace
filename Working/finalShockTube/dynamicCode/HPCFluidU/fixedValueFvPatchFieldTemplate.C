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
    // SHA1 = 9b8b4e26a42272db7722b37117fea651ec7f91c3
    //
    // unique function name that can be checked if the correct library version
    // has been loaded
    void HPCFluidU_9b8b4e26a42272db7722b37117fea651ec7f91c3(bool load)
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
    HPCFluidUFixedValueFvPatchVectorField
);


const char* const HPCFluidUFixedValueFvPatchVectorField::SHA1sum =
    "9b8b4e26a42272db7722b37117fea651ec7f91c3";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

HPCFluidUFixedValueFvPatchVectorField::
HPCFluidUFixedValueFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF)
{
    if (false)
    {
        Info<<"construct HPCFluidU sha1: 9b8b4e26a42272db7722b37117fea651ec7f91c3"
            " from patch/DimensionedField\n";
    }
}


HPCFluidUFixedValueFvPatchVectorField::
HPCFluidUFixedValueFvPatchVectorField
(
    const HPCFluidUFixedValueFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper)
{
    if (false)
    {
        Info<<"construct HPCFluidU sha1: 9b8b4e26a42272db7722b37117fea651ec7f91c3"
            " from patch/DimensionedField/mapper\n";
    }
}


HPCFluidUFixedValueFvPatchVectorField::
HPCFluidUFixedValueFvPatchVectorField
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
        Info<<"construct HPCFluidU sha1: 9b8b4e26a42272db7722b37117fea651ec7f91c3"
            " from patch/dictionary\n";
    }
}


HPCFluidUFixedValueFvPatchVectorField::
HPCFluidUFixedValueFvPatchVectorField
(
    const HPCFluidUFixedValueFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf)
{
    if (false)
    {
        Info<<"construct HPCFluidU sha1: 9b8b4e26a42272db7722b37117fea651ec7f91c3"
            " as copy\n";
    }
}


HPCFluidUFixedValueFvPatchVectorField::
HPCFluidUFixedValueFvPatchVectorField
(
    const HPCFluidUFixedValueFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF)
{
    if (false)
    {
        Info<<"construct HPCFluidU sha1: 9b8b4e26a42272db7722b37117fea651ec7f91c3 "
            "as copy/DimensionedField\n";
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

HPCFluidUFixedValueFvPatchVectorField::
~HPCFluidUFixedValueFvPatchVectorField()
{
    if (false)
    {
        Info<<"destroy HPCFluidU sha1: 9b8b4e26a42272db7722b37117fea651ec7f91c3\n";
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void HPCFluidUFixedValueFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (false)
    {
        Info<<"updateCoeffs HPCFluidU sha1: 9b8b4e26a42272db7722b37117fea651ec7f91c3\n";
    }

//{{{ begin code
    #line 188 "/home/anthonygay1812/OpenFOAM/Working/finalShockTube/system/fluid/codeDict.HPCFluidU"
const fvMesh& fluidMesh = db().parent().lookupObject<fvMesh>("fluid");
    const fvMesh& partMesh = db().parent().lookupObject<fvMesh>("particle");

    const volVectorField& UFluid = db().lookupObject<volVectorField>("U");     
    const volVectorField& UFluidOld=UFluid.oldTime();
   
    const volVectorField& momentum = partMesh.lookupObject<volVectorField>("momentum");
    const volScalarField& rhoM = partMesh.lookupObject<volScalarField>("rhoM");

    volVectorField U
            (
                IOobject
                (
                    "U",
                    partMesh.time().timeName(),
                    partMesh,
                    IOobject::NO_READ
                ),
                momentum/rhoM
            );

    const label fluidIdx=fluidMesh.nCells()-1;
    const label particleIdx=0;
    const scalar w=0.4; 
    
    Info<< nl << "FLUID VEL BC";
    Info<< nl << "Particle velocity: " << U[particleIdx];
        Info<< nl << "Fluid velocity : " << UFluid[fluidIdx];

    operator==((1-w)*UFluidOld[fluidIdx]+w*(0.5*UFluid[fluidIdx]+0.5*U[particleIdx]));
//}}} end code

    this->fixedValueFvPatchField<vector>::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

