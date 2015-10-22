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

#include "addedMassModel.H"


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::addedMassModel::update()
{
    // Set all time dependent variables for previous time step
    timePrev_ = time_.value();
    alphaPrev_ = alpha_;
    normalRelVelPrev_ = normalRelVel_;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::addedMassModel::addedMassModel
(
    const Time& time,
    scalar chordLength
)
:
    time_(time),
    chordLength_(chordLength),
    timePrev_(time.value()),
    nNewTimes_(0),
    alpha_(0.0),
    alphaPrev_(0.0),
    chordwiseRelVel_(0.0),
    normalRelVel_(0.0),
    normalRelVelPrev_(0.0)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::addedMassModel>
Foam::addedMassModel::New()
{
    return autoPtr<addedMassModel>(new addedMassModel);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::addedMassModel::~addedMassModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::addedMassModel::correct
(
    scalar& liftCoefficient,
    scalar& dragCoefficient,
    scalar& momentCoefficient,
    scalar alphaRad,
    scalar chordwiseRelVel,
    scalar normalRelVel
)
{
    scalar pi = Foam::constant::mathematical::pi;
    scalar time = time_.value();
    scalar deltaT = time_.deltaT().value();
    chordwiseRelVel_ = chordwiseRelVel;
    normalRelVel_ = normalRelVel;
    scalar relVelMagSqr = magSqr(normalRelVel_) + magSqr(chordwiseRelVel_);
    
    // Update previous values if time has changed
    if (time != timePrev_)
    {
        nNewTimes_++;
        if (nNewTimes_ > 1) 
        {
            update();
        }
    }
    
    if (nNewTimes_ <= 1)
    {
        alpha_ = alphaRad;
        alphaPrev_ = alpha_;
    }
    
    // Calculate added mass coefficients for a flat plate
    scalar alphaDot = (alpha_ - alphaPrev_)/deltaT;
    scalar ct = pi/8.0*chordLength_*alphaDot*normalRelVel_/relVelMagSqr;
    scalar normVelDot = (normalRelVel_ - normalRelVelPrev_)/deltaT;
    scalar cn = pi/8.0*chordLength_*normVelDot/relVelMagSqr;
    // Moment coefficient is at quarter chord
    scalar cm = cn/4.0 - pi/8.0*normalRelVel_*chordwiseRelVel_/relVelMagSqr;
    
    // Modify lift and drag coefficients
    liftCoefficient += cn*cos(alpha_) + ct*sin(alpha_);
    dragCoefficient += cn*sin(alpha_) - ct*cos(alpha_);
    momentCoefficient += cm;
}


// ************************************************************************* //
