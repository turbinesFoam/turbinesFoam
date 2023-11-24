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

#include "LeishmanBeddoes3G.H"
#include "addToRunTimeSelectionTable.H"
#include "unitConversion.H"

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
    if (mag(alphaEquiv_) > 2*Foam::constant::mathematical::pi)
    {
        alphaEquiv_ =  fmod(alphaEquiv_, 2*Foam::constant::mathematical::pi);
    }
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

    // Calculate the impulsive moment coefficient
    lambdaM_ = 3*pi/16*(alpha_ + c_/(4*magU_)*deltaAlpha_/deltaT_)
             + pi/16*c_/magU_*deltaAlpha_/deltaT_;
    J_ = JPrev_*exp(-deltaT_/TI_)
       + (lambdaM_ - lambdaMPrev_)*exp(-deltaT_/(2.0*TI_));
    CMI_ = -4.0/M_*J_;

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

    if (debug)
    {
        Info<< "    lambdaL: " << lambdaL_ << endl;
        Info<< "    lambdaLPrev: " << lambdaLPrev_ << endl;
        Info<< "    TI: " << TI_ << endl;
        Info<< "    H: " << H_ << endl;
    }
}


void Foam::fv::LeishmanBeddoes3G::calcS1S2
(
    scalar B,
    scalar C,
    scalar D
)
{
    LeishmanBeddoes::calcS1S2(B, C, D);
}


void Foam::fv::LeishmanBeddoes3G::calcSeparated()
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

    // Evaluate vortex tracking time
    if (not stalledPrev_)
    {
        tau_ = 0.0;
    }
    else
    {
        if (tau_ == tauPrev_)
        {
            tau_ = tauPrev_ + deltaS_;
        }
    }

    // Modify Tf time constant if necessary
    scalar Tf = Tf_;
    if (tau_ > Tvl_)
    {
        Tf = 0.5*Tf_;
    }

    // Calculate dynamic separation point
    scalar pi = Foam::constant::mathematical::pi;
    DF_ = DFPrev_*exp(-deltaS_/Tf)
        + (fPrime_ - fPrimePrev_)*exp(-deltaS_/(2*Tf));
    fDoublePrime_ = fPrime_ - DF_;

    if (tau_ >= 0 and tau_ <= Tvl_)
    {
        Vx_ = pow((sin(pi*tau_/(2.0*Tvl_))), 1.5);
    }
    else if (tau_ > Tvl_)
    {
        Vx_ = pow((cos(pi*(tau_ - Tvl_)/Tv_)), 2);
    }
    if (mag(alpha_) < mag(alphaPrev_))
    {
        Vx_ = 0.0;
    }

    // Calculate the separation point and limit to [0, 1]
    f3G_ = fDoublePrime_ - DF_*Vx_;
    if (f3G_ < 0)
    {
        f3G_ = 0.0;
    }
    else if (f3G_ > 1)
    {
        f3G_ = 1.0;
    }

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

    // Total normal force coefficient does not have CNV contribution since
    // this is included in the Vx term
    CN_ = CNF_;

    // Calculate moment coefficient (needs to be verified in reference)
    scalar m = cmFitExponent_;
    scalar cmf = (K0_ + K1_*(1 - fDoublePrime_)
               + K2_*sin(pi*Foam::pow(fDoublePrime_, m)))*CNC_
               + profileData_.zeroLiftMomentCoeff();
    scalar cmv = 0.2*(1.0 - cos(pi*tau_/Tvl_))*CNV_;
    CM_ = cmf + cmv + CMI_;
}


void Foam::fv::LeishmanBeddoes3G::update()
{
    LeishmanBeddoes::update();
    ZPrev_ = Z_;
    etaLPrev_ = etaL_;
    HPrev_ = H_;
    lambdaLPrev_ = lambdaL_;
    JPrev_ = J_;
    lambdaMPrev_ = lambdaM_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::LeishmanBeddoes3G::LeishmanBeddoes3G
(
    const dictionary& dict,
    const word& modelName,
    const Time& time,
    profileData& profileData
)
:
    LeishmanBeddoes(dict, modelName, time, profileData),
    Z_(0.0),
    ZPrev_(0.0),
    etaL_(0.0),
    etaLPrev_(0.0),
    A3_(coeffs_.lookupOrDefault("A3", 0.5)),
    T1_(coeffs_.lookupOrDefault("T1", 20.0)),
    T2_(coeffs_.lookupOrDefault("T2", 4.5)),
    H_(0.0),
    HPrev_(0.0),
    lambdaL_(0.0),
    lambdaLPrev_(0.0),
    J_(0.0),
    JPrev_(0.0),
    lambdaM_(0.0),
    lambdaMPrev_(0.0),
    CMC_(0.0),
    CMI_(0.0)
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
