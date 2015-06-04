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

#include "LeishmanBeddoesSGC.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(LeishmanBeddoesSGC, 0);
    addToRunTimeSelectionTable
    (
        dynamicStallModel, 
        LeishmanBeddoesSGC,
        dictionary
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fv::LeishmanBeddoesSGC::calcAlphaEquiv()
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


void Foam::fv::LeishmanBeddoesSGC::calcUnsteady()
{
    LeishmanBeddoes3G::calcUnsteady();
    
    // Calculate lagged angle of attack
    scalar DAlpha = alphaPrimePrev_*exp(-deltaS_/TAlpha_) 
                  + (alpha_ - alphaPrev_)*exp(-deltaS_/(2.0*TAlpha_));
    alphaPrime_ = alpha_ - DAlpha;
    
    // Calculate reduced pitch rate
    r_ = deltaAlpha_/deltaT_*c_/(2.0*magU_);
    
    // Calculate alphaDS0
    scalar pi = Foam::constant::mathematical::pi;
    scalar dAlphaDS = alphaDS0DiffDeg_/180.0*pi;
    alphaDS0_ = alphaSS_ + dAlphaDS;
    
    if (r_ >= r0_)
    {
        alphaCrit_ = alphaDS0_;
    }
    else
    {
        alphaCrit_ = alphaSS_ + (alphaDS0_ - alphaSS_)*r_/r0_;
    }
    
    stalled_ = (mag(alphaPrime_) > alphaCrit_);
}


void Foam::fv::LeishmanBeddoesSGC::calcSeparated()
{
    // Calculate trailing-edge separation point
    if (mag(alphaPrime_) < alpha1_)
    {
        fPrime_ = 1.0 - 0.4*exp((mag(alphaPrime_) - alpha1_)/S1_);
    }
    else
    {
        fPrime_ = 0.02 + 0.58*exp((alpha1_ - mag(alphaPrime_))/S2_);
    }
    
    // Modify Tf time constant if necessary
    scalar Tf = Tf_;
    if (tau_ > Tvl_) Tf = 0.5*Tf_;
    
    // Calculate dynamic separation point
    scalar pi = Foam::constant::mathematical::pi;
    DF_ = DFPrev_*exp(-deltaS_/Tf) 
        + (fPrime_ - fPrimePrev_)*exp(-deltaS_/(2*Tf));
    fDoublePrime_ = mag(fPrime_ - DF_);
    if (tau_ > 0 and tau_ <= Tvl_)
    {
        Vx_ = pow((sin(pi*tau_/(2.0*Tvl_))), 1.5);
    }
    else if (tau_ > Tvl_)
    {
        Vx_ = pow((cos(pi*(tau_ - Tvl_)/Tv_)), 2);
    }
    if (mag(alpha_) < mag(alphaPrev_)) Vx_ = 0.0;
    f3G_ = mag(fDoublePrime_ - DF_*Vx_);
    
    // Calculate normal force coefficient including dynamic separation point
    CNF_ = CNAlpha_*alphaEquiv_*pow(((1.0 + sqrt(f3G_))/2.0), 2) 
         + CNI_;
    
    // Calculate tangential force coefficient
    if (fDoublePrime_ < fCrit_)
    {
        CT_ = eta_*CNAlpha_*alphaEquiv_*alphaEquiv_*pow(fDoublePrime_, 1.5);
    }
    else
    {
        CT_ = eta_*CNAlpha_*alphaEquiv_*alphaEquiv_*sqrt(fDoublePrime_);
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
            CV_ = CNC_*(1.0 - pow(((1.0 + sqrt(fDoublePrime_))/2.0), 2));
            CNV_ = CNVPrev_*exp(-deltaS_/Tv) 
                 + (CV_ - CVPrev_)*exp(-deltaS_/(2.0*Tv));
        }
        else
        {
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


void Foam::fv::LeishmanBeddoesSGC::update()
{
    LeishmanBeddoes3G::update();
    alphaPrimePrev_ = alphaPrime_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::LeishmanBeddoesSGC::LeishmanBeddoesSGC
(
    const dictionary& dict,
    const word& modelName,
    const Time& time
)
:
    LeishmanBeddoes3G(dict, modelName, time),
    TAlpha_(coeffs_.lookupOrDefault("TAlpha", 6.30)),
    r0_(coeffs_.lookupOrDefault("r0", 0.01)),
    B1_(coeffs_.lookupOrDefault("B1", 0.5)),
    E0_(coeffs_.lookupOrDefault("E0", 0.15)),
    alphaDS0DiffDeg_(coeffs_.lookupOrDefault("alphaDS0DiffDeg", 3.6))
{
    Tv_ = coeffs_.lookupOrDefault("Tv", 11.0);
    Tvl_ = coeffs_.lookupOrDefault("Tvl", 9.0);
    eta_ = coeffs_.lookupOrDefault("eta", 0.975);
    
    if (debug)
    {
        Info<< modelName << " dynamic stall model created" << endl
            << "    Coeffs:" << endl << coeffs_ << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::LeishmanBeddoesSGC::~LeishmanBeddoesSGC()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


// ************************************************************************* //
