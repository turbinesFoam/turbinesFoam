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

#include "dynamicStallModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const dataType Foam::dynamicStallModel::staticData();


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dynamicStallModel::dynamicStallModel()
:
    baseClassName(),
    data_()
{}


Foam::dynamicStallModel::dynamicStallModel(const dataType& data)
:
    baseClassName(),
    data_(data)
{}


Foam::dynamicStallModel::dynamicStallModel(const dynamicStallModel&)
:
    baseClassName(),
    data_()
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::dynamicStallModel>
Foam::dynamicStallModel::New()
{
    return autoPtr<dynamicStallModel>(new dynamicStallModel);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dynamicStallModel::~dynamicStallModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

void Foam::dynamicStallModel::operator=(const dynamicStallModel& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn("Foam::dynamicStallModel::operator=(const Foam::dynamicStallModel&)")
            << "Attempted assignment to self"
            << abort(FatalError);
    }
}

// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //


// ************************************************************************* //
