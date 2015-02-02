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
    Foam::LeishmanBeddoes

Description
    Leishman-Beddoes dynamic stall model for use with actuatorLineElement.

SourceFiles
    LeishmanBeddoes.C

\*---------------------------------------------------------------------------*/

#ifndef LeishmanBeddoes_H
#define LeishmanBeddoes_H

#include "dynamicStallModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class LeishmanBeddoes Declaration
\*---------------------------------------------------------------------------*/

class LeishmanBeddoes
:
    public dynamicStallModel
{
    // Private data

    // Private Member Functions

        //- Disallow default bitwise copy construct
        LeishmanBeddoes(const LeishmanBeddoes&);

        //- Disallow default bitwise assignment
        void operator=(const LeishmanBeddoes&);


public:

    // Static data members

        //- Static data staticData
        static const dataType staticData;


    // Constructors

        //- Construct null
        LeishmanBeddoes();

        //- Construct from components
        LeishmanBeddoes(const dataType& data);

        //- Construct from Istream
        LeishmanBeddoes(Istream&);

        //- Construct as copy
        LeishmanBeddoes(const LeishmanBeddoes&);


    // Selectors

        //- Select null constructed
        static autoPtr<LeishmanBeddoes> New();


    //- Destructor
    ~LeishmanBeddoes();


    // Member Functions

        // Access

        // Check

        // Edit

        // Write


    // Member Operators

        void operator=(const LeishmanBeddoes&);


    // Friend Functions

    // Friend Operators

    // IOstream Operators

        friend Istream& operator>>(Istream&, LeishmanBeddoes&);
        friend Ostream& operator<<(Ostream&, const LeishmanBeddoes&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "LeishmanBeddoesI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
