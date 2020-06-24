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
#line 179 "/home/anthonygay1812/OpenFOAM/Working/finalShockTube/system/fluid/codeDict.HPCFluidT"
#include "constants.H"
  using namespace Foam::constant;
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
    // SHA1 = 2dc3145a9cfa6beb53a9a17d24491db4278714fa
    //
    // unique function name that can be checked if the correct library version
    // has been loaded
    void HPCFluidT_2dc3145a9cfa6beb53a9a17d24491db4278714fa(bool load)
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
    fvPatchScalarField,
    HPCFluidTFixedValueFvPatchScalarField
);


const char* const HPCFluidTFixedValueFvPatchScalarField::SHA1sum =
    "2dc3145a9cfa6beb53a9a17d24491db4278714fa";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

HPCFluidTFixedValueFvPatchScalarField::
HPCFluidTFixedValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF)
{
    if (false)
    {
        Info<<"construct HPCFluidT sha1: 2dc3145a9cfa6beb53a9a17d24491db4278714fa"
            " from patch/DimensionedField\n";
    }
}


HPCFluidTFixedValueFvPatchScalarField::
HPCFluidTFixedValueFvPatchScalarField
(
    const HPCFluidTFixedValueFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper)
{
    if (false)
    {
        Info<<"construct HPCFluidT sha1: 2dc3145a9cfa6beb53a9a17d24491db4278714fa"
            " from patch/DimensionedField/mapper\n";
    }
}


HPCFluidTFixedValueFvPatchScalarField::
HPCFluidTFixedValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<scalar>(p, iF, dict)
{
    if (false)
    {
        Info<<"construct HPCFluidT sha1: 2dc3145a9cfa6beb53a9a17d24491db4278714fa"
            " from patch/dictionary\n";
    }
}


HPCFluidTFixedValueFvPatchScalarField::
HPCFluidTFixedValueFvPatchScalarField
(
    const HPCFluidTFixedValueFvPatchScalarField& ptf
)
:
    fixedValueFvPatchField<scalar>(ptf)
{
    if (false)
    {
        Info<<"construct HPCFluidT sha1: 2dc3145a9cfa6beb53a9a17d24491db4278714fa"
            " as copy\n";
    }
}


HPCFluidTFixedValueFvPatchScalarField::
HPCFluidTFixedValueFvPatchScalarField
(
    const HPCFluidTFixedValueFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(ptf, iF)
{
    if (false)
    {
        Info<<"construct HPCFluidT sha1: 2dc3145a9cfa6beb53a9a17d24491db4278714fa "
            "as copy/DimensionedField\n";
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

HPCFluidTFixedValueFvPatchScalarField::
~HPCFluidTFixedValueFvPatchScalarField()
{
    if (false)
    {
        Info<<"destroy HPCFluidT sha1: 2dc3145a9cfa6beb53a9a17d24491db4278714fa\n";
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void HPCFluidTFixedValueFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (false)
    {
        Info<<"updateCoeffs HPCFluidT sha1: 2dc3145a9cfa6beb53a9a17d24491db4278714fa\n";
    }

//{{{ begin code
    #line 114 "/home/anthonygay1812/OpenFOAM/Working/finalShockTube/system/fluid/codeDict.HPCFluidT"
const fvMesh& fluidMesh = db().parent().lookupObject<fvMesh>("fluid");
    const fvMesh& partMesh = db().parent().lookupObject<fvMesh>("particle");

    const volScalarField& TFluid = db().lookupObject<volScalarField>("T");     
    const volScalarField& TFluidOld=TFluid.oldTime();

    const volVectorField& momentum = partMesh.lookupObject<volVectorField>("momentum");
    const volScalarField& rhoN = partMesh.lookupObject<volScalarField>("rhoN");
    const volScalarField& rhoM = partMesh.lookupObject<volScalarField>("rhoM");
    const volScalarField& linearKE = partMesh.lookupObject<volScalarField>("linearKE");
    const volScalarField& internalE = partMesh.lookupObject<volScalarField>("internalE");
    const volScalarField& iDof = partMesh.lookupObject<volScalarField>("iDof");

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
            
volScalarField translationalT
    (
        IOobject
        (
            "TTranslationalPart",
            partMesh.time().timeName(),
            partMesh,
            IOobject::NO_READ
        ),
        2.0/(3.0*physicoChemical::k.value()*rhoN)*(linearKE - 0.5*rhoM*(U & U))
    );

    volScalarField overallT
            (
                IOobject
                (
                    "overallT",
                    partMesh.time().timeName(),
                    partMesh,
                    IOobject::NO_READ
                ),
                2.0/(physicoChemical::k.value()*(3.0*rhoN + iDof))
              *(linearKE - 0.5*rhoM*(U & U) + internalE)
            );


    const label fluidIdx=fluidMesh.nCells()-1;
    const label particleIdx=0;
    const scalar w=0.4; 
    
     Info<< nl << "FLUID TEMP BC";
    Info<< nl << "Particle overallT: " << overallT[particleIdx];
    Info<< nl << "Particle TRans T: " << translationalT[particleIdx];
    Info<< nl << "Fluid Temp: " << TFluid[fluidIdx] << nl;

    operator==((1-w)*TFluidOld[fluidIdx]+w*(0.5*TFluid[fluidIdx]+0.5*translationalT[particleIdx]));
//}}} end code

    this->fixedValueFvPatchField<scalar>::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

