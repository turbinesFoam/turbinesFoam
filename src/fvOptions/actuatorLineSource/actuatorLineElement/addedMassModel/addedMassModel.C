/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author(s)
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of turbinesFoam, which is based on OpenFOAM.

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
    scalar chordLength,
    const label& debug
)
:
    time_(time),
    chordLength_(chordLength),
    debug(debug),
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
Foam::addedMassModel::New
(
    const Time& time,
    scalar chordLength,
    const label& debug
)
{
    return autoPtr<addedMassModel>
    (
        new addedMassModel
        (
            time,
            chordLength,
            debug
        )
    );
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
    alpha_ = alphaRad;
    scalar alphaDot = (alpha_ - alphaPrev_)/deltaT;
    scalar normVelDot = (normalRelVel_ - normalRelVelPrev_)/deltaT;

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
    scalar cc = pi/8.0*chordLength_*alphaDot*normalRelVel_/relVelMagSqr;
    scalar cn = -pi/8.0*chordLength_*normVelDot/relVelMagSqr;
    // Moment coefficient is at quarter chord
    scalar cm = -cn/4.0 - pi/8.0*normalRelVel_*chordwiseRelVel_/relVelMagSqr;

    if (debug)
    {
        Info<< endl << "Added mass quantities:" << endl;
        Info<< "    alphaDot: " << alphaDot << endl;
        Info<< "    normVelDot: " << normVelDot << endl;
        Info<< "    cn: " << cn << endl;
        Info<< "    cc: " << cc << endl;
        Info<< "    cm: " << cm << endl << endl;
        Info<< "Coefficients before added mass correction:" << endl;
        Info<< "    cl: " << liftCoefficient  << endl;
        Info<< "    cd: " << dragCoefficient << endl;
        Info<< "    cm: " << momentCoefficient << endl << endl;
    }

    // Modify lift and drag coefficients
    liftCoefficient += cn*cos(alpha_) + cc*sin(alpha_);
    dragCoefficient += cn*sin(alpha_) - cc*cos(alpha_);
    momentCoefficient += cm;

    if (debug)
    {
        Info<< "Coefficients after added mass correction:" << endl;
        Info<< "    cl: " << liftCoefficient << endl;
        Info<< "    cd: " << dragCoefficient << endl;
        Info<< "    cm: " << momentCoefficient << endl << endl;
    }
}


// ************************************************************************* //
