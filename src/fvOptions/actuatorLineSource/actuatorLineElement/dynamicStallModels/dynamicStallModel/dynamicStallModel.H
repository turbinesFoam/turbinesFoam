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
    Foam::dynamicStallModel

Description

SourceFiles
    dynamicStallModelI.H
    dynamicStallModel.C
    dynamicStallModelIO.C

\*---------------------------------------------------------------------------*/

#ifndef dynamicStallModel_H
#define dynamicStallModel_H

#include ".H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;

// Forward declaration of friend functions and operators
class dynamicStallModel;
Istream& operator>>(Istream&, dynamicStallModel&);
Ostream& operator<<(Ostream&, const dynamicStallModel&);


/*---------------------------------------------------------------------------*\
                         Class dynamicStallModel Declaration
\*---------------------------------------------------------------------------*/

class dynamicStallModel
:
    public baseClassName
{
    // Private data

        //- Description of data_
        dataType data_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        dynamicStallModel(const dynamicStallModel&);

        //- Disallow default bitwise assignment
        void operator=(const dynamicStallModel&);


public:

    // Static data members

        //- Static data staticData
        static const dataType staticData;


    // Constructors

        //- Construct null
        dynamicStallModel();

        //- Construct from components
        dynamicStallModel(const dataType& data);

        //- Construct from Istream
        dynamicStallModel(Istream&);

        //- Construct as copy
        dynamicStallModel(const dynamicStallModel&);


    // Selectors

        //- Select null constructed
        static autoPtr<dynamicStallModel> New();


    //- Destructor
    ~dynamicStallModel();


    // Member Functions

        // Access

        // Check

        // Edit

        // Write


    // Member Operators

        void operator=(const dynamicStallModel&);


    // Friend Functions

    // Friend Operators

    // IOstream Operators

        friend Istream& operator>>(Istream&, dynamicStallModel&);
        friend Ostream& operator<<(Ostream&, const dynamicStallModel&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "dynamicStallModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
