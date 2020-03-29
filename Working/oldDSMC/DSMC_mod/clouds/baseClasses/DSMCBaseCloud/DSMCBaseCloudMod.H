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
    Foam::DSMCBaseCloudMod

Description
    Virtual abstract base class for templated DSMCCloud

SourceFiles
    DSMCBaseCloudMod.C

\*---------------------------------------------------------------------------*/

#ifndef DSMCBaseCloudMod_H
#define DSMCBaseCloudMod_H

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class DSMCBaseCloudMod Declaration
\*---------------------------------------------------------------------------*/

class DSMCBaseCloudMod
{
public:

    //- Runtime type information
    TypeName("DSMCBaseCloudMod");

    // Constructors

        //- Null constructor
        DSMCBaseCloudMod();

        //- Disallow default bitwise copy construction
        DSMCBaseCloudMod(const DSMCBaseCloudMod&) = delete;


    //- Destructor
    virtual ~DSMCBaseCloudMod();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const DSMCBaseCloudMod&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
