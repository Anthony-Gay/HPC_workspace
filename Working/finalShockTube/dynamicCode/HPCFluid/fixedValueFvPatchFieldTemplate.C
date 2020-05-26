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
#line 82 "/home/anthonygay1812/OpenFOAM/Working/finalShockTube/system/fluid/codeDict.HPCFluid"
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
    // SHA1 = 120a66b9aef436830054c2b9d97943264fd6a35a
    //
    // unique function name that can be checked if the correct library version
    // has been loaded
    void HPCFluid_120a66b9aef436830054c2b9d97943264fd6a35a(bool load)
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
    HPCFluidFixedValueFvPatchScalarField
);


const char* const HPCFluidFixedValueFvPatchScalarField::SHA1sum =
    "120a66b9aef436830054c2b9d97943264fd6a35a";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

HPCFluidFixedValueFvPatchScalarField::
HPCFluidFixedValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF)
{
    if (false)
    {
        Info<<"construct HPCFluid sha1: 120a66b9aef436830054c2b9d97943264fd6a35a"
            " from patch/DimensionedField\n";
    }
}


HPCFluidFixedValueFvPatchScalarField::
HPCFluidFixedValueFvPatchScalarField
(
    const HPCFluidFixedValueFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper)
{
    if (false)
    {
        Info<<"construct HPCFluid sha1: 120a66b9aef436830054c2b9d97943264fd6a35a"
            " from patch/DimensionedField/mapper\n";
    }
}


HPCFluidFixedValueFvPatchScalarField::
HPCFluidFixedValueFvPatchScalarField
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
        Info<<"construct HPCFluid sha1: 120a66b9aef436830054c2b9d97943264fd6a35a"
            " from patch/dictionary\n";
    }
}


HPCFluidFixedValueFvPatchScalarField::
HPCFluidFixedValueFvPatchScalarField
(
    const HPCFluidFixedValueFvPatchScalarField& ptf
)
:
    fixedValueFvPatchField<scalar>(ptf)
{
    if (false)
    {
        Info<<"construct HPCFluid sha1: 120a66b9aef436830054c2b9d97943264fd6a35a"
            " as copy\n";
    }
}


HPCFluidFixedValueFvPatchScalarField::
HPCFluidFixedValueFvPatchScalarField
(
    const HPCFluidFixedValueFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(ptf, iF)
{
    if (false)
    {
        Info<<"construct HPCFluid sha1: 120a66b9aef436830054c2b9d97943264fd6a35a "
            "as copy/DimensionedField\n";
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

HPCFluidFixedValueFvPatchScalarField::
~HPCFluidFixedValueFvPatchScalarField()
{
    if (false)
    {
        Info<<"destroy HPCFluid sha1: 120a66b9aef436830054c2b9d97943264fd6a35a\n";
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void HPCFluidFixedValueFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (false)
    {
        Info<<"updateCoeffs HPCFluid sha1: 120a66b9aef436830054c2b9d97943264fd6a35a\n";
    }

//{{{ begin code
    #line 20 "/home/anthonygay1812/OpenFOAM/Working/finalShockTube/system/fluid/codeDict.HPCFluid"
const fvMesh& fluidMeshRef = db().parent().lookupObject<fvMesh>("fluid");
    const fvMesh& partMeshRef = db().parent().lookupObject<fvMesh>("particle");

    const volScalarField& pFluid = db().lookupObject<volScalarField>("p");     
    const volScalarField& pFluidOld=pFluid.oldTime();

    const volVectorField& momentumRef = partMeshRef.lookupObject<volVectorField>("momentum");
    const volScalarField& rhoNRef = partMeshRef.lookupObject<volScalarField>("rhoN");
    const volScalarField& rhoMRef = partMeshRef.lookupObject<volScalarField>("rhoM");
    const volScalarField& linearKERef = partMeshRef.lookupObject<volScalarField>("linearKE");

  volVectorField URef
    (
        IOobject
        (
            "UPart",
            partMeshRef.time().timeName(),
            partMeshRef,
            IOobject::NO_READ
        ),
        momentumRef/rhoMRef
    );

  volScalarField translationalTRef
    (
        IOobject
        (
            "TTranslationalPart",
            partMeshRef.time().timeName(),
            partMeshRef,
            IOobject::NO_READ
        ),
        2.0/(3.0*physicoChemical::k.value()*rhoNRef)*(linearKERef - 0.5*rhoMRef*(URef & URef))
    );

    volScalarField pParticle
    (
        IOobject
        (
            "pParticle",
            partMeshRef.time().timeName(),
            partMeshRef,
            IOobject::NO_READ
        ),
        physicoChemical::k.value()*rhoNRef*translationalTRef
    );
    
    const label fluidIdx=fluidMeshRef.nCells()-1;;
    const label particleIdx=0;
    const scalar w=0.2; 
    
    Info<< nl << "FLUID BC";
    Info<< nl << "Expected Pressure: 1.13E2";
    Info<< nl << "Calculated Pressure: "<< ((1-w)*pFluidOld[fluidIdx]+w*(0.5*pFluid[fluidIdx]+0.5*pParticle[particleIdx]));
    Info<< nl << "Difference: "<< 1.13E2-((1-w)*pFluidOld[fluidIdx]+w*(0.5*pFluid[fluidIdx]+0.5*pParticle[particleIdx])) << endl;

    operator==((1-w)*pFluidOld[fluidIdx]+w*(0.5*pFluid[fluidIdx]+0.5*pParticle[particleIdx]));
//}}} end code

    this->fixedValueFvPatchField<scalar>::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

