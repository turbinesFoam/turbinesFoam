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

#include "LeishmanBeddoesSD.H"
#include "addToRunTimeSelectionTable.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(LeishmanBeddoesSD, 0);
    addToRunTimeSelectionTable
    (
        dynamicStallModel,
        LeishmanBeddoesSD,
        dictionary
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fv::LeishmanBeddoesSD::calcUnsteady()
{
    if (not alphaAttachedCorrection_)
    {
        alphaEquiv_ = alpha_;
    }

    LeishmanBeddoes3G::calcUnsteady();

    if (not alphaAttachedCorrection_)
    {
        CNI_ = 0;
    }

    // First part of the downwind size, accelerate stall for cross flow turbines
    scalar TAlpha = TAlpha_;
    if (crossFlowTurbine_)
    {
        if (alphaEquiv_ < 0 && deltaAlpha_ < 0)
        {
            TAlpha *= 0.1;
        }
    }

    // Calculate lagged angle of attack
    DAlpha_ = DAlphaPrev_*exp(-deltaSPrev_/TAlpha)
            + (alpha_ - alphaPrev_)*exp(-deltaS_/(2.0*TAlpha));
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


void Foam::fv::LeishmanBeddoesSD::calcSeparated()
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

    // First part of the downwind size, accelerate stall for cross flow turbines
    scalar Tf = Tf_;
    if (crossFlowTurbine_)
    {
        if (alphaEquiv_ < 0 && deltaAlpha_ < 0)
        {
            Tf *= 0.1;
        }
    }

    // Calculate dynamic separation point
    scalar pi = Foam::constant::mathematical::pi;
    DF_ = DFPrev_*exp(-deltaSPrev_/Tf)
        + (fPrime_ - fPrimePrev_)*exp(-deltaS_/(2*Tf));
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

    if (fDoublePrime_ < 0.0)
    {
        fDoublePrime_ = 0;
    }

    // Calculate normal force coefficient including dynamic separation point
    CNF_ = CNAlpha_*alphaEquiv_*pow(((1.0 + sqrt(fDoublePrime_))/2.0), 2)
         + CNI_;

    // Calculate tangential force coefficient
    CT_ = eta_*CNAlpha_*alphaEquiv_*alphaEquiv_*(sqrt(fDoublePrime_)
        - E0_*pow(fDoublePrime_, 1.0/Tv_));


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

    if (f < 0.0)
    {
        f = 0.0;
    }

    CTStatic_ = eta_*CNAlpha_*alphaEquiv_*alphaEquiv_*(sqrt(f)
              - E0_*pow(f, 1.0/Tv_));

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


void Foam::fv::LeishmanBeddoesSD::update()
{
    LeishmanBeddoesSGC::update();
    deltaSPrev_ = deltaS_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::LeishmanBeddoesSD::LeishmanBeddoesSD
(
    const dictionary& dict,
    const word& modelName,
    const Time& time,
    profileData& profileData
)
:
    LeishmanBeddoesSGC(dict, modelName, time, profileData),
    CTCorrection_(coeffs_.lookupOrDefault("CTCorrection", true)),
    alphaAttachedCorrection_
    (
        coeffs_.lookupOrDefault("alphaAttachedCorrection", true)
    ),
    crossFlowTurbine_(coeffs_.lookupOrDefault("crossFlowTurbine", false)),
    deltaSPrev_(0.0),
    CTStatic_(0.0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::LeishmanBeddoesSD::~LeishmanBeddoesSD()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::fv::LeishmanBeddoesSD::correct
(
    scalar magU,
    scalar alphaDeg,
    scalar& cl,
    scalar& cd,
    scalar& cm
)
{
    LeishmanBeddoesSGC::correct
    (
        magU,
        alphaDeg,
        cl,
        cd,
        cm
    );

    // CTCorrection is a correction to ensure that the model reduces to static
    // values when pitch rate and angle of attack is low enough.
    if (CTCorrection_)
    {
        const scalar CC_corr_const = 1;
        const scalar scaleStart = CC_corr_const * 0.5 * alphaSS_;
        const scalar scaleEnd = CC_corr_const * alphaSS_;
        const scalar rScaleStart = 0.01;
        const scalar rScaleEnd = 0.02;
        scalar alphaDeg = radToDeg(alpha_);

        scalar scaleFactorAlpha =
            (scaleEnd - fabs(alphaDeg))/(scaleEnd-scaleStart);
        if (scaleFactorAlpha < 0)
        {
            scaleFactorAlpha = 0;
        }
        if (scaleFactorAlpha > 1)
        {
            scaleFactorAlpha = 1;
        }

        scalar scaleFactor_r =
            (rScaleEnd - abs(r_))/(rScaleEnd-rScaleStart);
        if (scaleFactor_r < 0)
        {
            scaleFactor_r = 0;
        }
        if (scaleFactor_r > 1)
        {
            scaleFactor_r = 1;
        }
        scalar scaleFactor = scaleFactorAlpha*scaleFactor_r;

        scalar CTProfileData = profileData_.chordwiseCoefficient(alphaDeg);
        scalar CTVal = CT_ + (CTProfileData - CTStatic_)*scaleFactor;

        cl = CN_*cos(alpha_) + CTVal*sin(alpha_);
        cd = CN_*sin(alpha_) - CTVal*cos(alpha_) + CD0_*(1 - scaleFactor);
    }
}

// ************************************************************************* //
