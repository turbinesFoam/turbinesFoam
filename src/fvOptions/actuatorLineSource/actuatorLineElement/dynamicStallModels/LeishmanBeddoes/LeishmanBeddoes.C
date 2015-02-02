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

\*---------------------------------------------------------------------------*/

#include "LeishmanBeddoes.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const dataType Foam::LeishmanBeddoes::staticData();


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::LeishmanBeddoes::LeishmanBeddoes()
:
    baseClassName(),
    data_()
{}


Foam::LeishmanBeddoes::LeishmanBeddoes(const dataType& data)
:
    baseClassName(),
    data_(data)
{}


Foam::LeishmanBeddoes::LeishmanBeddoes(const LeishmanBeddoes&)
:
    baseClassName(),
    data_()
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::LeishmanBeddoes>
Foam::LeishmanBeddoes::New()
{
    return autoPtr<LeishmanBeddoes>(new LeishmanBeddoes);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::LeishmanBeddoes::~LeishmanBeddoes()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

void Foam::LeishmanBeddoes::operator=(const LeishmanBeddoes& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn("Foam::LeishmanBeddoes::operator=(const Foam::LeishmanBeddoes&)")
            << "Attempted assignment to self"
            << abort(FatalError);
    }
}

// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //


// ************************************************************************* //
