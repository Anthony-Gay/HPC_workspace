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
    // SHA1 = 63ac19f166ab850d826a1946c60e72e3aaae4679
    //
    // unique function name that can be checked if the correct library version
    // has been loaded
    void HPCFluid_63ac19f166ab850d826a1946c60e72e3aaae4679(bool load)
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
    "63ac19f166ab850d826a1946c60e72e3aaae4679";


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
        Info<<"construct HPCFluid sha1: 63ac19f166ab850d826a1946c60e72e3aaae4679"
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
        Info<<"construct HPCFluid sha1: 63ac19f166ab850d826a1946c60e72e3aaae4679"
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
        Info<<"construct HPCFluid sha1: 63ac19f166ab850d826a1946c60e72e3aaae4679"
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
        Info<<"construct HPCFluid sha1: 63ac19f166ab850d826a1946c60e72e3aaae4679"
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
        Info<<"construct HPCFluid sha1: 63ac19f166ab850d826a1946c60e72e3aaae4679 "
            "as copy/DimensionedField\n";
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

HPCFluidFixedValueFvPatchScalarField::
~HPCFluidFixedValueFvPatchScalarField()
{
    if (false)
    {
        Info<<"destroy HPCFluid sha1: 63ac19f166ab850d826a1946c60e72e3aaae4679\n";
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
        Info<<"updateCoeffs HPCFluid sha1: 63ac19f166ab850d826a1946c60e72e3aaae4679\n";
    }

//{{{ begin code
    #line 30 "/home/anthonygay1812/OpenFOAM/tutorials/tutorials/compressible/rhoCentralFoam/shockTube/0/p.boundaryField.sides"
operator==(0.5*sin(3.14*this->db().time().value())/2);
//}}} end code

    this->fixedValueFvPatchField<scalar>::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

