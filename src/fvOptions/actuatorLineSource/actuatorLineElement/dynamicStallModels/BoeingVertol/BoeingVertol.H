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
    Foam::BoeingVertol

Description

SourceFiles
    BoeingVertolI.H
    BoeingVertol.C
    BoeingVertolIO.C

\*---------------------------------------------------------------------------*/

#ifndef BoeingVertol_H
#define BoeingVertol_H

#include ".H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;

// Forward declaration of friend functions and operators
class BoeingVertol;
Istream& operator>>(Istream&, BoeingVertol&);
Ostream& operator<<(Ostream&, const BoeingVertol&);


/*---------------------------------------------------------------------------*\
                         Class BoeingVertol Declaration
\*---------------------------------------------------------------------------*/

class BoeingVertol
:
    public baseClassName
{
    // Private data

        //- Description of data_
        dataType data_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        BoeingVertol(const BoeingVertol&);

        //- Disallow default bitwise assignment
        void operator=(const BoeingVertol&);


public:

    // Static data members

        //- Static data staticData
        static const dataType staticData;


    // Constructors

        //- Construct null
        BoeingVertol();

        //- Construct from components
        BoeingVertol(const dataType& data);

        //- Construct from Istream
        BoeingVertol(Istream&);

        //- Construct as copy
        BoeingVertol(const BoeingVertol&);


    // Selectors

        //- Select null constructed
        static autoPtr<BoeingVertol> New();


    //- Destructor
    ~BoeingVertol();


    // Member Functions

        // Access

        // Check

        // Edit

        // Write


    // Member Operators

        void operator=(const BoeingVertol&);


    // Friend Functions

    // Friend Operators

    // IOstream Operators

        friend Istream& operator>>(Istream&, BoeingVertol&);
        friend Ostream& operator<<(Ostream&, const BoeingVertol&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "BoeingVertolI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
