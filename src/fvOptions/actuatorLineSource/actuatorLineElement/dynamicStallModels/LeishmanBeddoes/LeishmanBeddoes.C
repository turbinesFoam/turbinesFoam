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
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(LeishmanBeddoes, 0);
    addToRunTimeSelectionTable
    (
        dynamicStallModel, 
        LeishmanBeddoes,
        dictionary
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fv::LeishmanBeddoes::update()
{
    timePrev_ = time_;
    alphaPrev_ = alpha_;
    XPrev_ = X_;
    YPrev_ = Y_;
    deltaAlphaPrev_ = alphaPrev_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::LeishmanBeddoes::LeishmanBeddoes
(
    const dictionary& dict,
    const word& modelName,
    const scalar startTime
)
:
    dynamicStallModel(dict, modelName, startTime),
    
    A1_(coeffs_.lookupOrDefault("A1", 0.3)),
    A2_(coeffs_.lookupOrDefault("A2", 0.7)),
    b1_(coeffs_.lookupOrDefault("b1", 0.14)),
    b2_(coeffs_.lookupOrDefault("b2", 0.53)),
    a_(coeffs_.lookupOrDefault("speedOfSound", 1e12))
{
    dict_.lookup("chordLength") >> c_;
    
    if (debug)
    {
        Info<< modelName << " dynamic stall model created" << endl
            << "    Coeffs:" << endl << coeffs_ << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::LeishmanBeddoes::~LeishmanBeddoes()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::fv::LeishmanBeddoes::correct
(
    scalar alphaDeg, 
    scalar& cl, 
    scalar& cd
)
{
}

void Foam::fv::LeishmanBeddoes::correct
(
    scalar time,
    scalar magU,
    scalar alphaDeg,
    scalar& cl,
    scalar& cd,
    List<scalar> alphaDegList,
    List<scalar> clList,
    List<scalar> cdList
)
{
    time_ = time;
    alpha_ = alphaDeg/180*Foam::constant::mathematical::pi;
    
    if (time_ != timePrev_)
    {
        update();
    }
}

// ************************************************************************* //
