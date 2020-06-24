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
#line 104 "/home/anthonygay1812/OpenFOAM/Working/finalShockTube/system/fluid/codeDict.HPCFluidP"
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
    // SHA1 = a6d23aedfcb9d8349035819c70b7760d0e9994dd
    //
    // unique function name that can be checked if the correct library version
    // has been loaded
    void HPCFluidP_a6d23aedfcb9d8349035819c70b7760d0e9994dd(bool load)
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
    HPCFluidPFixedValueFvPatchScalarField
);


const char* const HPCFluidPFixedValueFvPatchScalarField::SHA1sum =
    "a6d23aedfcb9d8349035819c70b7760d0e9994dd";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

HPCFluidPFixedValueFvPatchScalarField::
HPCFluidPFixedValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF)
{
    if (false)
    {
        Info<<"construct HPCFluidP sha1: a6d23aedfcb9d8349035819c70b7760d0e9994dd"
            " from patch/DimensionedField\n";
    }
}


HPCFluidPFixedValueFvPatchScalarField::
HPCFluidPFixedValueFvPatchScalarField
(
    const HPCFluidPFixedValueFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper)
{
    if (false)
    {
        Info<<"construct HPCFluidP sha1: a6d23aedfcb9d8349035819c70b7760d0e9994dd"
            " from patch/DimensionedField/mapper\n";
    }
}


HPCFluidPFixedValueFvPatchScalarField::
HPCFluidPFixedValueFvPatchScalarField
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
        Info<<"construct HPCFluidP sha1: a6d23aedfcb9d8349035819c70b7760d0e9994dd"
            " from patch/dictionary\n";
    }
}


HPCFluidPFixedValueFvPatchScalarField::
HPCFluidPFixedValueFvPatchScalarField
(
    const HPCFluidPFixedValueFvPatchScalarField& ptf
)
:
    fixedValueFvPatchField<scalar>(ptf)
{
    if (false)
    {
        Info<<"construct HPCFluidP sha1: a6d23aedfcb9d8349035819c70b7760d0e9994dd"
            " as copy\n";
    }
}


HPCFluidPFixedValueFvPatchScalarField::
HPCFluidPFixedValueFvPatchScalarField
(
    const HPCFluidPFixedValueFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(ptf, iF)
{
    if (false)
    {
        Info<<"construct HPCFluidP sha1: a6d23aedfcb9d8349035819c70b7760d0e9994dd "
            "as copy/DimensionedField\n";
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

HPCFluidPFixedValueFvPatchScalarField::
~HPCFluidPFixedValueFvPatchScalarField()
{
    if (false)
    {
        Info<<"destroy HPCFluidP sha1: a6d23aedfcb9d8349035819c70b7760d0e9994dd\n";
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void HPCFluidPFixedValueFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (false)
    {
        Info<<"updateCoeffs HPCFluidP sha1: a6d23aedfcb9d8349035819c70b7760d0e9994dd\n";
    }

//{{{ begin code
    #line 20 "/home/anthonygay1812/OpenFOAM/Working/finalShockTube/system/fluid/codeDict.HPCFluidP"
const fvMesh& fluidMesh = db().parent().lookupObject<fvMesh>("fluid");
    const fvMesh& partMesh = db().parent().lookupObject<fvMesh>("particle");

    const volScalarField& pFluid = db().lookupObject<volScalarField>("p");     
    const volScalarField& pFluidOld=pFluid.oldTime();

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
            "UPart",
            partMesh.time().timeName(),
            partMesh,
            IOobject::NO_READ
        ),
        momentum/rhoM
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

    volScalarField pParticle
    (
        IOobject
        (
            "pParticle",
            partMesh.time().timeName(),
            partMesh,
            IOobject::NO_READ
        ),
        physicoChemical::k.value()*rhoN*translationalT //https://courses.lumenlearning.com/physics/chapter/13-3-the-ideal-gas-law/
    );
    
    const label fluidIdx=fluidMesh.nCells()-1;
    const label particleIdx=0;
    const scalar w=0.4; 
    
    Info<< nl << "FLUID PRESSURE BC";
    Info<< nl << "Particle pressure formula 1: " << pParticle[particleIdx];
    Info<< nl << "Particle pressure formula 2: " << physicoChemical::k.value()*rhoN[particleIdx]*translationalT[particleIdx];
    Info<< nl << "Particle pressure from rhoM: " << rhoM[particleIdx]*overallT[particleIdx]*296.8;
        Info<< nl << "Particle pressure from rhoM T_trans: " << rhoM[particleIdx]*translationalT[particleIdx]*296.8;
    Info<< nl << "Particle overallT: " << overallT[particleIdx];
    Info<< nl << "Particle TRans T: " << translationalT[particleIdx] << nl;
  

    operator==((1-w)*pFluidOld[fluidIdx]+w*(0.5*pFluid[fluidIdx]+0.5*pParticle[particleIdx]));
//}}} end code

    this->fixedValueFvPatchField<scalar>::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

