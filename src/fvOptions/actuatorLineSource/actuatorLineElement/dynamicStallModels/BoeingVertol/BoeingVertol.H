/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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
    Foam::fv::BoeingVertol

Description
    Boeing-Vertol dynamic stall model for use with actuatorLineElement.

SourceFiles
    BoeingVertol.C

\*---------------------------------------------------------------------------*/

#ifndef BoeingVertol_H
#define BoeingVertol_H

#include "dynamicStallModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                         Class BoeingVertol Declaration
\*---------------------------------------------------------------------------*/

class BoeingVertol
:
    public dynamicStallModel
{


public:

    //- Runtime type information
    TypeName("BoeingVertol");

    // Constructors

        //- Construct null
        BoeingVertol
        (
            const dictionary& dict,
            const word& modelName
        );

        //- Construct as copy
        BoeingVertol(const BoeingVertol&);


    // Selectors

        //- Select null constructed
        static autoPtr<BoeingVertol> New
        (
            const dictionary& dict,
            const word& modelName
        );


    //- Destructor
    ~BoeingVertol();


    // Member Functions

        // Access

        // Check

        // Edit

        // Write

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
