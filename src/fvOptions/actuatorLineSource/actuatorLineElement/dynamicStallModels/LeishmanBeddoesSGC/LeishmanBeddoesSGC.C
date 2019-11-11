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

#include "LeishmanBeddoesSGC.H"
#include "addToRunTimeSelectionTable.H"
#include "unitConversion.H"

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

void Foam::fv::LeishmanBeddoesSGC::calcUnsteady()
{
    LeishmanBeddoes3G::calcUnsteady();

    // Calculate lagged angle of attack
    DAlpha_ = DAlphaPrev_*exp(-deltaS_/TAlpha_)
            + (alpha_ - alphaPrev_)*exp(-deltaS_/(2.0*TAlpha_));
    alphaPrime_ = alpha_ - DAlpha_;

    // Calculate reduced pitch rate
    r_ = deltaAlpha_/deltaT_*c_/(2.0*magU_);

    // Calculate alphaDS0
    scalar pi = Foam::constant::mathematical::pi;
    scalar dAlphaDS = alphaDS0DiffDeg_/180.0*pi;
    alphaDS0_ = alphaSS_ + dAlphaDS;

    if (mag(r_) >= r0_)
    {
        alphaCrit_ = alphaDS0_;
    }
    else
    {
        alphaCrit_ = alphaSS_ + (alphaDS0_ - alphaSS_)*mag(r_)/r0_;
    }

    stalled_ = (mag(alphaPrime_) > alphaCrit_);

    if (debug)
    {
        Info<< "    Reduced pitch rate: " << r_ << endl;
        Info<< "    alphaCrit: " << alphaCrit_ << endl;
    }
}


void Foam::fv::LeishmanBeddoesSGC::calcSeparated()
{
    // Calculate lagged trailing-edge separation point
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

    // Calculate dynamic separation point
    scalar pi = Foam::constant::mathematical::pi;
    DF_ = DFPrev_*exp(-deltaS_/Tf_)
        + (fPrime_ - fPrimePrev_)*exp(-deltaS_/(2*Tf_));
    fDoublePrime_ = fPrime_ - DF_;

    // Calculate vortex modulation parameter
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

    // Calculate normal force coefficient including dynamic separation point
    CNF_ = CNAlpha_*alphaEquiv_*pow(((1.0 + sqrt(fDoublePrime_))/2.0), 2)
         + CNI_;

    // Calculate tangential force coefficient
    CT_ = eta_*CNAlpha_*alphaEquiv_*alphaEquiv_*(sqrt(fDoublePrime_) - E0_);

    // Calculate static trailing-edge separation point
    scalar f;
    if (mag(alpha_) < alpha1_)
    {
        f = 1.0 - 0.4*exp((mag(alpha_) - alpha1_)/S1_);
    }
    else
    {
        f = 0.02 + 0.58*exp((alpha1_ - mag(alpha_))/S2_);
    }

    // Evaluate vortex lift contributions
    CNV_ = B1_*(fDoublePrime_ - f)*Vx_;

    // Total normal force coefficient is the combination of that from
    // circulatory effects, impulsive effects, dynamic separation, and vortex
    // lift
    CN_ = CNF_ + CNV_;

    // Calculate moment coefficient
    scalar m = cmFitExponent_;
    scalar cmf = (K0_ + K1_*(1 - fDoublePrime_)
               + K2_*sin(pi*Foam::pow(fDoublePrime_, m)))*CNC_
               + profileData_.zeroLiftMomentCoeff();
    scalar cmv = B2_*(1.0 - cos(pi*tau_/Tvl_))*CNV_;
    CM_ = cmf + cmv + CMI_;

    if (debug)
    {
        Info<< "    Vx: " << Vx_ << endl;
        Info<< "    f: " << f << endl;
        Info<< "    cmf: " << cmf << endl;
        Info<< "    cmv: " << cmv << endl;
    }
}


void Foam::fv::LeishmanBeddoesSGC::update()
{
    LeishmanBeddoes3G::update();
    alphaPrimePrev_ = alphaPrime_;
    DAlphaPrev_ = DAlpha_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::LeishmanBeddoesSGC::LeishmanBeddoesSGC
(
    const dictionary& dict,
    const word& modelName,
    const Time& time,
    profileData& profileData
)
:
    LeishmanBeddoes3G(dict, modelName, time, profileData),
    TAlpha_(coeffs_.lookupOrDefault("TAlpha", 6.30)),
    alphaPrimePrev_(0.0),
    alphaCrit_(17.0),
    r_(0.0),
    r0_(coeffs_.lookupOrDefault("r0", 0.01)),
    B1_(coeffs_.lookupOrDefault("B1", 0.5)),
    B2_(coeffs_.lookupOrDefault("B2", 0.2)),
    E0_(coeffs_.lookupOrDefault("E0", 0.15)),
    alphaDS0DiffDeg_(coeffs_.lookupOrDefault("alphaDS0DiffDeg", 3.6)),
    DAlpha_(0.0),
    DAlphaPrev_(0.0)
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
