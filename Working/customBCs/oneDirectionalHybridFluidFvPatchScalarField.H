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

Class
    Foam::oneDirectionalHybridFluidFvPatchScalarField

Description
    This boundary condition provides a total temperature condition.

Usage
    \table
        Property     | Description             | Required    | Default value
        U            | Velocity field name     | no          | U
        phi          | Flux field name         | no          | phi
        psi          | Compressibility field name | no       | thermo:psi
        gamma        | ratio of specific heats (Cp/Cv) | yes |
        T0           | reference temperature   | yes         |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            oneDirectionalHybridFluid;
        T0              uniform 300;
    }
    \endverbatim

SourceFiles
    oneDirectionalHybridFluidFvPatchScalarField.C

See also
    Foam::fixedValueFvPatchField

\*---------------------------------------------------------------------------*/

#ifndef oneDirectionalHybridFluidFvPatchScalarField_H
#define oneDirectionalHybridFluidFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class oneDirectionalHybridFluidFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class oneDirectionalHybridFluidFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private Data
         label fluidIdx_;
         label particleIdx_;
         scalar w_;

public:

    //- Runtime type information
    TypeName("oneDirectionalHybridFluid");


    // Constructors

        //- Construct from patch and internal field
        oneDirectionalHybridFluidFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        oneDirectionalHybridFluidFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given oneDirectionalHybridFluidFvPatchScalarField
        //  onto a new patch
        oneDirectionalHybridFluidFvPatchScalarField
        (
            const oneDirectionalHybridFluidFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy constructor
        oneDirectionalHybridFluidFvPatchScalarField
        (
            const oneDirectionalHybridFluidFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new oneDirectionalHybridFluidFvPatchScalarField(*this)
            );
        }

        //- Copy constructor setting internal field reference
        oneDirectionalHybridFluidFvPatchScalarField
        (
            const oneDirectionalHybridFluidFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new oneDirectionalHybridFluidFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        // Access

            //- Return relaxation coeff
            scalar relaxationCoeff() const
            {
                return w_;
            }

            label particleIndex() const
            {
                return particleIdx_;
            }; 
         
            label fluidIndex() const
            {
                return fluidIdx_;
            }; 

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            //  Used to update fields following mesh topology change
            virtual void autoMap(const fvPatchFieldMapper&);

            //- Reverse map the given fvPatchField onto this fvPatchField
            //  Used to reconstruct fields
            virtual void rmap(const fvPatchScalarField&, const labelList&);


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
