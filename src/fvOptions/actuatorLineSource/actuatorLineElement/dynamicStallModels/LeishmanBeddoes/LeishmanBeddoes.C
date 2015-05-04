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

void Foam::fv::LeishmanBeddoes::evalStaticData
(
    scalar alphaDeg,
    List<scalar> alphaDegList,
    List<scalar> clList,
    List<scalar> cdList
)
{
    // Create lists for normal and chordwise coefficients
    scalar pi = Foam::constant::mathematical::pi;
    scalar alphaRad = alphaDeg/180*pi;
    List<scalar> alphaRadList = alphaDegList/180*pi;
    List<scalar> cnList = clList*sin(alphaRadList) - cdList*cos(alphaRadList);
    
    // Calculate lift slope CNAlpha
    scalar cn0 = interpolate(0, alphaDegList, cnList);
    scalar cn5 = interpolate(5, alphaDegList, cnList);
    CNAlpha_ = (cn5 - cn0)/(5/180*pi);
    
    // Calculate critical normal force coefficient CN1, where the slope of the
    // drag coefficient curve slope first breaks 0.02 per degree
    forAll(alphaDegList, i)
    {
        scalar alpha = alphaDegList[i];
        scalar cd0, cd1;
        if (alpha > 4 && alpha < 25)
        {
            cd1 = interpolate(alpha, alphaDegList, cdList);
            cd0 = interpolate(alpha + 0.5, alphaDegList, cdList);
            if ((cd1 - cd0)/0.5 > 0.02)
            {
                CN1_ = cnList[i];
                break;
            }
        }
    }
    
    // Calculate alpha1
    scalar f = 0.7;
    alpha1_ = CN1_/CNAlpha_/pow((1 + sqrt(f))/2, 2);
    
    // Calculate S1 or S2, depending on whether alpha is above or below alpha1
    if (abs(alphaRad) < alpha1_)
    {
        S1_ = (abs(alphaRad) - alpha1_)/log((f - 1)/(-0.3));
    }
    else if (abs(alphaRad) >= alpha1_)
    {
        S2_ = (alpha1_ - abs(alphaRad))/log((f - 0.04)/(0.66));
    }
    
    // Calculate CD0
    CD0_ = interpolate(0, alphaDegList, cdList);
}


void Foam::fv::LeishmanBeddoes::update()
{
    timePrev_ = time_;
    alphaPrev_ = alpha_;
    XPrev_ = X_;
    YPrev_ = Y_;
    deltaAlphaPrev_ = deltaAlpha_;
    DPrev_ = D_;
    DPPrev_ = DP_;
    CNPPrev_ = CNP_;
    DFPrev_ = DF_;
    fPrimePrev_ = fPrime_;
    CVPrev_ = CV_;
    CNVPrev_ = CNV_;
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
    time_ = startTime;
    timePrev_ = startTime;
    
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
    M_ = magU/a_;
    deltaAlpha_ = alpha_ - alphaPrev_;
    
    // Only calculate deltaT if time has changed
    if (time != timePrev_)
    {
        deltaT_ = time_ - timePrev_;
    }
    
    deltaS_ = 2*magU*deltaT_/c_;
    
    if (debug)
    {
        Info<< "Leishman-Beddoes dynamic stall model correcting" << endl;
        Info<< "deltaT: " << deltaT_ << endl;
        Info<< "deltaAlpha: " << deltaAlpha_ << endl;
    }
    
    if (time_ != timePrev_)
    {
        update();
    }
}

// ************************************************************************* //
