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
    List<scalar> alphaRadList(alphaDegList.size());
    List<scalar> cnList(clList.size());
    
    forAll(alphaDegList, i)
    {
        alphaRadList[i] = alphaDegList[i]/180*pi;
        cnList[i] = clList[i]*sin(alphaRadList[i]) 
                  - cdList[i]*cos(alphaRadList[i]);
    }
    
    // Calculate lift slope CNAlpha
    scalar cn0 = interpolate(0.0, alphaDegList, cnList);
    scalar cn5 = interpolate(5.0, alphaDegList, cnList);
    scalar dAlpha = 5.0/180.0*pi;
    CNAlpha_ = (cn5 - cn0)/dAlpha;
    
    // Calculate critical normal force coefficient CN1, where the slope of the
    // drag coefficient curve slope first breaks 0.02 per degree
    forAll(alphaDegList, i)
    {
        scalar alpha = alphaDegList[i];
        scalar cd0, cd1, slope, dAlpha;
        if (alpha > 4 && alpha < 25)
        {
            cd1 = cdList[i + 1];
            cd0 = cdList[i];
            dAlpha = alphaDegList[i + 1] - alpha;
            slope = (cd1 - cd0)/dAlpha;
            if (slope > 0.02)
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


void Foam::fv::LeishmanBeddoes::calcUnsteady()
{
    // Calculate the equivalent angle of attack
    scalar beta = 1 - M_;
    X_ = XPrev_*exp(-b1_*beta*deltaS_) 
       + A1_*deltaAlpha_*exp(b1_*beta*deltaS_/2);
    Y_ = YPrev_*exp(-b2_*beta*deltaS_) 
       + A2_*deltaAlpha_*exp(b2_*beta*deltaS_/2);
    alphaEquiv_ = alpha_ - X_ - Y_;
    
    // Calculate the circulatory normal force coefficient
    CNC_ = CNAlpha_*alphaEquiv_;
    
    // Calculate the impulsive normal force coefficient
    scalar pi = Foam::constant::mathematical::pi;
    scalar kAlpha = 0.75/(1 - M_ + pi*(1 - M_*M_)*M_*M_*(A1_*b1_ + A2_*b2_));
    TI_ = c_/a_;
    D_ = DPrev_*exp(-deltaT_/(kAlpha*TI_)) 
       - ((deltaAlpha_ - deltaAlphaPrev_)/deltaT_)
       *exp(-deltaT_/(2*kAlpha*TI_));
    CNI_ = 4*kAlpha*TI_/M_*(deltaAlpha_/deltaT_ - D_);
    
    // Calculate total normal force coefficient
    CNP_ = CNC_ + CNI_;
    
    // Apply first-order lag to normal force coefficient
    DP_ = DPPrev_*exp(-deltaS_/Tp_) + (CNP_ - CNPPrev_)*exp(-deltaS_/(2*Tp_));
    CNPrime_ = CNP_ - DP_;
    
    // Calculate lagged angle of attack
    alphaPrime_ = CNPrime_/CNAlpha_;
    
    // Set stalled switch
    stalled_ = (CNPrime_ > CN1_);
}


void Foam::fv::LeishmanBeddoes::calcSeparated()
{
    // Calculate trailing-edge separation point
    if (abs(alphaPrime_) < alpha1_)
    {
        fPrime_ = 1.0 - 0.3*exp((abs(alphaPrime_) - alpha1_)/S1_);
    }
    else if (abs(alphaPrime_) >= alpha1_)
    {
        fPrime_ = 0.04 + 0.66*exp((alpha1_ - abs(alphaPrime_))/S2_);
    }
    
    // Calculate dynamic separation point
    DF_ = DFPrev_*exp(-deltaS_/Tf_) 
        + (fPrime_ - fPrimePrev_)*exp(-deltaS_/(2*Tf_));
    fDoublePrime_ = fPrime_ - DF_;
    
    // Calculate normal force coefficient for dynamic separation point
    CNF_ = CNAlpha_*alphaEquiv_*pow(((1 + sqrt(fDoublePrime_))/2), 2) + CNI_;
    
    // Calculate tangential force coefficient
    if (fDoublePrime_ < 0.7)
    {
        CT_ = eta_*CNAlpha_*alphaEquiv_*alphaEquiv_*sqrt(fDoublePrime_);
    }
    else
    {
        CT_ = eta_*CNAlpha_*alphaEquiv_*alphaEquiv_*pow(fDoublePrime_, 1.5);
    }
    
    // Compute vortex shedding process if stalled
    if (stalled_)
    {
        // Evaluate vortex tracking time
        if (not stalledPrev_) tau_ = 0.0;
        else tau_ = tauPrev_ + deltaS_;
        
        // Evaluate vortex lift contributions
        if (tau_ < Tvl_)
        {
            CV_ = CNC_*(1 - pow(((1 + sqrt(fDoublePrime_))/2), 2));
            CNV_ = CNVPrev_*exp(-deltaS_/Tv_) 
                 + (CV_ - CVPrev_)*exp(-deltaS_/(2*Tv_));
        }
        else
        {
            CNV_ = CNVPrev_*exp(-deltaS_/Tv_);
        }
    }
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
    stalledPrev_ = stalled_;
    tauPrev_ = tau_;
    firstTime_ = false;
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
    a_(coeffs_.lookupOrDefault("speedOfSound", 1e12)),
    eta_(coeffs_.lookupOrDefault("eta", 0.95)),
    Tp_(coeffs_.lookupOrDefault("Tp", 1.7)),
    Tf_(coeffs_.lookupOrDefault("Tf", 3.0)),
    Tv_(coeffs_.lookupOrDefault("Tv", 6.0)),
    Tvl_(coeffs_.lookupOrDefault("Tvl", 7.0))
{
    dict_.lookup("chordLength") >> c_;
    time_ = startTime;
    timePrev_ = startTime;
    firstTime_ = true;
    
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
    
    if (firstTime_) deltaAlpha_ = 0.0;
    else deltaAlpha_ = alpha_ - alphaPrev_;
    
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
    
    evalStaticData(alphaDeg, alphaDegList, clList, cdList);
    calcUnsteady();
    calcSeparated();
    
    if (stalled_)
    {
        CN_ = CNF_ + CNV_;
    }
    else
    {
        CN_ = CNP_;
    }
    
    // Modify lift and drag coefficients based on new normal force coefficient
    cl = CN_*cos(alpha_) + CT_*sin(alpha_);
    cd = CN_*sin(alpha_) - CT_*cos(alpha_) + CD0_;
    
    if (time_ != timePrev_)
    {
        update();
    }
}

// ************************************************************************* //
