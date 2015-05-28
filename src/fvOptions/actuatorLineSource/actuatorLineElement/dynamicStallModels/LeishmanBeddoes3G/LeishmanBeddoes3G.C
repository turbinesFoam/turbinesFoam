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

#include "LeishmanBeddoes3G.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(LeishmanBeddoes3G, 0);
    addToRunTimeSelectionTable
    (
        dynamicStallModel, 
        LeishmanBeddoes3G,
        dictionary
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fv::LeishmanBeddoes3G::calcAlphaEquiv()
{
    T3_ = 1.25*M_;
    scalar beta = 1 - M_*M_;
    X_ = XPrev_*exp(-beta*deltaS_/T1_) 
       + A1_*(etaL_ - etaLPrev_)*exp(-beta*deltaS_/(2.0*T1_));
    Y_ = YPrev_*exp(-beta*deltaS_/T2_) 
       + A2_*(etaL_ - etaLPrev_)*exp(-beta*deltaS_/(2.0*T2_));
    Z_ = ZPrev_*exp(-beta*deltaS_/T3_) 
       + A1_*(etaL_ - etaLPrev_)*exp(-beta*deltaS_/(2.0*T3_));
    etaL_ = alpha_ + c_/(2.0*magU_)*deltaAlpha_/deltaT_;
    alphaEquiv_ = etaL_ - X_ - Y_ - Z_;
}


void Foam::fv::LeishmanBeddoes3G::calcUnsteady()
{
    // Calculate the circulatory normal force coefficient
    CNC_ = CNAlpha_*alphaEquiv_;
    
    // Calculate the impulsive normal force coefficient
    scalar pi = Foam::constant::mathematical::pi;
    lambdaL_ = (pi/4.0)*(alpha_ + c_/(4.0*magU_)*deltaAlpha_/deltaT_);
    TI_ = c_/a_*(1.0 + 3.0*M_)/4.0;
    H_ = HPrev_*exp(-deltaT_/TI_) 
       + (lambdaL_ - lambdaLPrev_)*exp(-deltaT_/(2.0*TI_));
    CNI_ = 4.0/M_*H_;
    
    // Calculate total normal force coefficient
    CNP_ = CNC_ + CNI_;
    
    // Apply first-order lag to normal force coefficient
    DP_ = DPPrev_*exp(-deltaS_/Tp_) 
        + (CNP_ - CNPPrev_)*exp(-deltaS_/(2.0*Tp_));
    CNPrime_ = CNP_ - DP_;
    
    // Calculate lagged angle of attack
    alphaPrime_ = CNPrime_/CNAlpha_;
    
    // Set stalled switch
    stalled_ = (mag(CNPrime_) > CN1_);
}


void Foam::fv::LeishmanBeddoes3G::calcS1S2
(
    List<scalar> alphaDegList,
    List<scalar> clList,
    List<scalar> cdList
)
{
    scalar pi = Foam::constant::mathematical::pi;
    scalar sumY = 0.0;
    scalar sumXYLnY = 0.0;
    scalar sumXY = 0.0;
    scalar sumYLnY = 0.0;
    scalar sumX2Y = 0.0;
    scalar alphaLowerLimit;
    scalar alphaUpperLimit;
    if (mag(alphaPrime_) < alpha1_)
    {
        alphaLowerLimit = 0.0;
        alphaUpperLimit = alpha1_;
    }
    else
    {
        alphaLowerLimit = alpha1_ - 1e-3;
        alphaUpperLimit = pi/6.0;
    }
    forAll(alphaDegList, i)
    {
        scalar alphaRad = alphaDegList[i]/180.0*pi;
        scalar cn = clList[i]*cos(alphaRad) - cdList[i]*sin(alphaRad);
        scalar f = 1.0;
        if (alphaRad > alphaLowerLimit and alphaRad < alphaUpperLimit)
        {
            f = pow((sqrt(mag(cn)/CNAlpha_/mag(alphaRad))
                    *2.0 - 1.0), 2);
            scalar x;
            scalar y;
            if (mag(alphaPrime_) < alpha1_) 
            {
                x = mag(alphaRad) - alpha1_;
                y = (f - 1)/(-0.4);
            }
            else 
            {
                x = alpha1_ - mag(alphaRad);
                y = (f - 0.02)/0.58;
            }
            if (f > 0 and f < 1 and y > 0)
            {
                sumY += y;
                sumXYLnY += x*y*log(y);
                sumXY += x*y;
                sumYLnY += y*log(y);
                sumX2Y += x*x*y;
            }
        }
    }
    scalar b = (sumY*sumXYLnY - sumXY*sumYLnY)/(sumY*sumX2Y - sumXY*sumXY);
    if (mag(alphaPrime_) < alpha1_)
    {
        S1_ = 1.0/b;
        S2_ = 0.0;
    }
    else
    {
        S1_ = 0.0;
        S2_ = 1.0/b;
    }
    
    if (debug)
    {
        Info<< "    S1: " << S1_ << endl;
        Info<< "    S2: " << S2_ << endl;
    }
}


void Foam::fv::LeishmanBeddoes3G::calcSeparated()
{
    // Calculate trailing-edge separation point
    if (mag(alphaPrime_) < alpha1_)
    {
        fPrime_ = 1.0 - 0.3*exp((mag(alphaPrime_) - alpha1_)/S1_);
    }
    else
    {
        fPrime_ = 0.04 + 0.66*exp((alpha1_ - mag(alphaPrime_))/S2_);
    }
    
    // Modify Tf time constant if necessary
    scalar Tf = Tf_;
    if (tau_ > 0 and tau_ <= Tvl_) Tf = 0.5*Tf_;
    else if (tau_ > Tvl_ and tau_ <= 2*Tvl_) Tf = 4*Tf_;
    if (mag(alpha_) < mag(alphaPrev_) and mag(CNPrime_) < CN1_)
    {
        Tf = 0.5*Tf_;
    }
    
    // Calculate dynamic separation point
    DF_ = DFPrev_*exp(-deltaS_/Tf) 
        + (fPrime_ - fPrimePrev_)*exp(-deltaS_/(2*Tf));
    fDoublePrime_ = mag(fPrime_ - DF_);
    
    // Calculate normal force coefficient including dynamic separation point
    CNF_ = CNAlpha_*alphaEquiv_*pow(((1.0 + sqrt(fDoublePrime_))/2.0), 2) 
         + CNI_;
    
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
    // Evaluate vortex tracking time
    if (not stalledPrev_) tau_ = 0.0;
    else 
    {
        if (tau_ == tauPrev_)
        {
            tau_ = tauPrev_ + deltaS_;
        }
    }
    
    // Calculate Strouhal number time constant and set tau to zero to 
    // allow multiple vortex shedding
    scalar Tst = 2.0*(1.0 - fDoublePrime_)/0.19;
    if (tau_ > (Tvl_ + Tst)) tau_ = 0.0;
    
    // Evaluate vortex lift contributions, which are only nonzero if angle
    // of attack increased in magnitude
    if (mag(alpha_) > mag(alphaPrev_) and mag(alpha_ - alphaPrev_) > 0.01)
    {
        scalar Tv = Tv_;
        if (tau_ < Tvl_)
        {
            // Halve Tv if dAlpha/dt changes sign
            if (sign(deltaAlpha_) != sign(deltaAlphaPrev_)) Tv = 0.5*Tv_;
            CV_ = CNC_*(1.0 - pow(((1.0 + sqrt(fDoublePrime_))/2.0), 2));
            CNV_ = CNVPrev_*exp(-deltaS_/Tv) 
                 + (CV_ - CVPrev_)*exp(-deltaS_/(2.0*Tv));
        }
        else
        {
            Tv = 0.5*Tv_;
            CNV_ = CNVPrev_*exp(-deltaS_/Tv);
        }
    }
    else
    {
        CNV_ = 0.0;
    }

    // Total normal force coefficient is the combination of that from
    // circulatory effects, impulsive effects, dynamic separation, and vortex 
    // lift
    CN_ = CNF_ + CNV_;
}


void Foam::fv::LeishmanBeddoes3G::update()
{
    LeishmanBeddoes::update();
    ZPrev_ = Z_;
    etaLPrev_ = etaL_;
    HPrev_ = H_;
    lambdaLPrev_ = lambdaL_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::LeishmanBeddoes3G::LeishmanBeddoes3G
(
    const dictionary& dict,
    const word& modelName,
    const Time& time
)
:
    LeishmanBeddoes(dict, modelName, time),
    Z_(0.0),
    etaL_(0.0),
    A3_(coeffs_.lookupOrDefault("A3", 0.5)),
    T1_(coeffs_.lookupOrDefault("T1", 20.0)),
    T2_(coeffs_.lookupOrDefault("T2", 4.5)),
    H_(0.0),
    lambdaL_(0.0)
{
    fCrit_ = 0.6;
    Tv_ = coeffs_.lookupOrDefault("Tv", 10.0);
    Tvl_ = coeffs_.lookupOrDefault("Tvl", 8.0);
    A1_ = coeffs_.lookupOrDefault("A1", 0.165);
    A2_ = coeffs_.lookupOrDefault("A2", 0.335);
    
    if (debug)
    {
        Info<< modelName << " dynamic stall model created" << endl
            << "    Coeffs:" << endl << coeffs_ << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::LeishmanBeddoes3G::~LeishmanBeddoes3G()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


// ************************************************************************* //
