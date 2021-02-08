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
    // SHA1 = 703974d06aec963ff4176132dc0ce075ca2d440e
    //
    // unique function name that can be checked if the correct library version
    // has been loaded
    void HPCFluidU_703974d06aec963ff4176132dc0ce075ca2d440e(bool load)
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
    "703974d06aec963ff4176132dc0ce075ca2d440e";


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
        Info<<"construct HPCFluidU sha1: 703974d06aec963ff4176132dc0ce075ca2d440e"
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
        Info<<"construct HPCFluidU sha1: 703974d06aec963ff4176132dc0ce075ca2d440e"
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
        Info<<"construct HPCFluidU sha1: 703974d06aec963ff4176132dc0ce075ca2d440e"
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
        Info<<"construct HPCFluidU sha1: 703974d06aec963ff4176132dc0ce075ca2d440e"
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
        Info<<"construct HPCFluidU sha1: 703974d06aec963ff4176132dc0ce075ca2d440e "
            "as copy/DimensionedField\n";
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

HPCFluidUFixedValueFvPatchVectorField::
~HPCFluidUFixedValueFvPatchVectorField()
{
    if (false)
    {
        Info<<"destroy HPCFluidU sha1: 703974d06aec963ff4176132dc0ce075ca2d440e\n";
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
        Info<<"updateCoeffs HPCFluidU sha1: 703974d06aec963ff4176132dc0ce075ca2d440e\n";
    }

//{{{ begin code
    #line 184 "/home/anthonygay1812/OpenFOAM/Working/finalShockTube/system/fluid/codeDict.HPCFluidU"
const fvMesh& fluidMesh = db().parent().lookupObject<fvMesh>("fluid");
    const fvMesh& partMesh = db().parent().lookupObject<fvMesh>("particle"); 
    const volVectorField& fluidU = fluidMesh.lookupObject<volVectorField>("U"); 
    label fluidBoundPatchId=fluidMesh.boundaryMesh().findPatchID("fluidBound");
    const vector UFluidBoundaryOld
    (
        fluidU.oldTime().boundaryField()[fluidBoundPatchId][0]
    );

    const volVectorField& partMomentum = partMesh.lookupObject<volVectorField>("momentum");
    const volScalarField& partRhoM = partMesh.lookupObject<volScalarField>("rhoM");
    
    volVectorField partU
            (
                IOobject
                (
                    "uPart",
                    partMesh.time().timeName(),
                    partMesh,
                    IOobject::NO_READ
                ),
              partMomentum/partRhoM
            );
    

    const label fluidIdx=fluidMesh.nCells()-1;
    const label particleIdx=0;
    const scalar w=0.4; 
    const vector partUCell=partU[particleIdx];
    const vector partUX(partUCell.component(vector::X),0,0);

    operator==((1-w)*UFluidBoundaryOld+w*(0.5*fluidU[fluidIdx]+0.5*partU[particleIdx]));
//}}} end code

    this->fixedValueFvPatchField<vector>::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

