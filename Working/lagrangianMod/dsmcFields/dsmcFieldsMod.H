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
    Foam::functionObjects::dsmcFieldsMod

Description
    Calculate intensive fields:
    - UMean
    - translationalT
    - internalT
    - overallT
    from averaged extensive fields from a DSMC calculation.

See also
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    dsmcFieldsMod.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_dsmcFieldsMod_H
#define functionObjects_dsmcFieldsMod_H

#include "fvMeshFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                         Class dsmcFieldsMod Declaration
\*---------------------------------------------------------------------------*/

class dsmcFieldsMod
:
    public fvMeshFunctionObject
{
public:

    //- Runtime type information
    TypeName("dsmcFieldsMod");


    // Constructors

        //- Construct from Time and dictionary
        dsmcFieldsMod
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        dsmcFieldsMod(const dsmcFieldsMod&) = delete;


    //- Destructor
    virtual ~dsmcFieldsMod();


    // Member Functions

        //- Read the dsmcFieldsMod data
        virtual bool read(const dictionary&);

        //- Do nothing
        virtual bool execute();

        //- Calculate and write the DSMC fields
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const dsmcFieldsMod&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
