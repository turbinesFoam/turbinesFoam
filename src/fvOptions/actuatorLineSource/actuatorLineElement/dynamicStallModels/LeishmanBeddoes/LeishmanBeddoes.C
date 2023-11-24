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

#include "LeishmanBeddoes.H"
#include "addToRunTimeSelectionTable.H"
#include "simpleMatrix.H"
#include "interpolateUtils.H"
#include "unitConversion.H"

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

Foam::List<scalar> Foam::fv::LeishmanBeddoes::cnToF
(
    List<scalar> cnList,
    List<scalar> alphaRadList,
    bool limit
)
{
    List<scalar> fList(cnList.size());
    scalar alpha0 = degToRad(profileData_.zeroLiftAngleOfAttack());

    forAll(cnList, i)
    {
        if (alphaRadList[i] == 0.0)
        {
            alphaRadList[i] = SMALL;
        }
        fList[i] = magSqr
        (
            2*Foam::sqrt(mag(cnList[i]) /
            (CNAlpha_ * mag(alphaRadList[i] - alpha0)))
            - 1.0
        );
        if (limit)
        {
            if (fList[i] >= 1)
            {
                fList[i] = 1.0 - SMALL;
            }
            if (fList[i] <= 0)
            {
                fList[i] = SMALL;
            }
        }
    }
    return fList;
}


void Foam::fv::LeishmanBeddoes::calcAlphaEquiv()
{
    scalar beta = 1.0 - M_*M_;
    X_ = XPrev_*exp(-b1_*beta*deltaS_)
       + A1_*deltaAlpha_*exp(b1_*beta*deltaS_/2.0);
    Y_ = YPrev_*exp(-b2_*beta*deltaS_)
       + A2_*deltaAlpha_*exp(b2_*beta*deltaS_/2.0);
    alphaEquiv_ = alpha_ - X_ - Y_;
}


void Foam::fv::LeishmanBeddoes::evalStaticData()
{
    // Get static stall angle in radians
    if (coeffs_.found("alphaSSList"))
    {
        alphaSS_ = interpolateStaticParam("alphaSS");
    }
    else
    {
        alphaSS_ = profileData_.staticStallAngleRad();
    }

    // Get normal coefficient slope CNAlpha
    if (coeffs_.found("CNAlphaList"))
    {
        CNAlpha_ = interpolateStaticParam("CNAlpha");
    }
    else
    {
        CNAlpha_ = profileData_.normalCoeffSlope();
    }

    // Calculate CN1 using normal coefficient slope and critical f value
    // alpha1 is 87% of the static stall angle
    scalar f = fCrit_;
    if (coeffs_.found("alpha1List"))
    {
        alpha1_ = interpolateStaticParam("alpha1");
    }
    else
    {
        alpha1_ = alphaSS_*0.87;
    }
    if (coeffs_.found("CN1List"))
    {
        CN1_ = interpolateStaticParam("CN1");
    }
    else
    {
        CN1_ = CNAlpha_*alpha1_*pow((1.0 + sqrt(f))/2.0, 2);
    }

    if (debug)
    {
        Info<< "    Evaluating static foil data" << endl;
        scalar cn = CNAlpha_*alpha_;
        Info<< "    Static stall angle (deg): " << radToDeg(alphaSS_) << endl;
        Info<< "    Critical normal force coefficient: " << CN1_ << endl;
        Info<< "    Normal coefficient slope: " << CNAlpha_ << endl;
        Info<< "    Normal coefficient from slope: " << cn << endl;
    }

    // Get CD0
    if (coeffs_.found("CD0List"))
    {
        CD0_ = interpolateStaticParam("CD0");
    }
    else
    {
        CD0_ = profileData_.zeroLiftDragCoeff();
    }
    if (debug)
    {
        Info<< "    Cd_0: " << CD0_ << endl;
        Info<< "    alpha1: " << alpha1_ << endl;
    }

    // Calculate S1 and S2 constants for the separation point curve
    if (coeffs_.found("S1List") and coeffs_.found("S2List"))
    {
        S1_ = interpolateStaticParam("S1");
        S2_ = interpolateStaticParam("S2");
    }
    else
    {
        calcS1S2();
    }

    // Calculate the K1 and K2 constants for the moment coefficient
    if (coeffs_.found("K1List") and coeffs_.found("K2List"))
    {
        K1_ = interpolateStaticParam("K1");
        K2_ = interpolateStaticParam("K2");
    }
    else
    {
        calcK1K2();
    }
}


scalar Foam::fv::LeishmanBeddoes::interpolateStaticParam(word paramName)
{
    scalar Re = profileData_.Re();
    List<scalar> ReList = coeffs_.lookup("ReList");
    label interpIndex = interpolateUtils::binarySearch(ReList, Re);
    scalar interpFraction = interpolateUtils::getPart(Re, ReList, interpIndex);
    List<scalar> paramList = coeffs_.lookup(paramName + "List");

    return interpolateUtils::interpolate1D
    (
        interpFraction,
        paramList,
        interpIndex
    );
}


void Foam::fv::LeishmanBeddoes::calcUnsteady()
{
    // Calculate the circulatory normal force coefficient
    CNC_ = CNAlpha_*alphaEquiv_;

    // Calculate the impulsive normal force coefficient
    scalar pi = Foam::constant::mathematical::pi;
    scalar kAlpha = 0.75
                  / (1.0 - M_ + pi*(1.0 - M_*M_)*M_*M_*(A1_*b1_ + A2_*b2_));
    TI_ = c_/a_;
    D_ = DPrev_*exp(-deltaT_/(kAlpha*TI_))
       + ((deltaAlpha_ - deltaAlphaPrev_)/deltaT_)
       *exp(-deltaT_/(2.0*kAlpha*TI_));
    CNI_ = 4.0*kAlpha*TI_/M_*(deltaAlpha_/deltaT_ - D_);

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
        Info<< "    CNP: " << CNP_ << endl;
        Info<< "    CNPrime: " << CNPrime_ << endl;
        Info<< "    TI: " << TI_ << endl;
        Info<< "    D: " << D_ << endl;
        Info<< "    KAlpha: " << kAlpha << endl;
    }
}


void Foam::fv::LeishmanBeddoes::calcS1S2
(
    scalar B,
    scalar C,
    scalar D
)
{
    // Get subset of angle of attack list in radians
    List<scalar> alphaList = degToRad
    (
        profileData_.angleOfAttackList
        (
            0.5,
            radToDeg(alpha1_)
        )
    );

    // Least squares fit to find S1
    List<scalar> cnList = profileData_.normalCoefficientList
    (
        0.5,
        radToDeg(alpha1_)
    );
    List<scalar> f = cnToF(cnList, alphaList, true);
    List<scalar> x = alphaList - alpha1_;
    List<scalar> y = Foam::mag((f - 1.0)/(-B));
    scalar b = Foam::sum(x*y*Foam::log(y))/Foam::sum(x*x*y);
    S1_ = 1/b;

    // Least squares fit to find S2
    alphaList = degToRad
    (
        profileData_.angleOfAttackList
        (
            radToDeg(alpha1_) + 0.001,
            25.0
        )
    );
    cnList = profileData_.normalCoefficientList
    (
        radToDeg(alpha1_) + 0.001,
        25.0
    );
    f = cnToF(cnList, alphaList, true);
    x = alpha1_ - alphaList;
    y = Foam::mag((f - C)/D);
    b = Foam::sum(x*y*Foam::log(y))/Foam::sum(x*x*y);
    S2_ = 1/b;

    if (debug)
    {
        Info<< "    S1: " << S1_ << endl;
        Info<< "    S2: " << S2_ << endl;
    }
}


void Foam::fv::LeishmanBeddoes::calcK1K2()
{
    scalar pi = Foam::constant::mathematical::pi;
    List<scalar> alpha = degToRad(profileData_.angleOfAttackList(0.5, 25));
    List<scalar> cn = profileData_.normalCoefficientList(0.5, 25);
    List<scalar> cm = profileData_.momentCoefficientList(0.5, 25);
    List<scalar> f = cnToF(cn, alpha);
    scalar m = cmFitExponent_;
    simpleMatrix<scalar> A(2);
    A[0][0] = Foam::sum(magSqr(1.0 - f));
    A[0][1] = Foam::sum(sin(pi*Foam::pow(f, m)));
    A[1][0] = Foam::sum(sin(pi*Foam::pow(f, m))*(1.0 - f));
    A[1][1] = Foam::sum(magSqr(sin(pi*Foam::pow(f, m))));
    A.source()[0] = Foam::sum(cm/cn*(1.0 - f)) - K0_*Foam::sum((1.0 - f));
    A.source()[1] = Foam::sum(cm/cn*sin(pi*Foam::pow(f, m)))
                  - K0_*Foam::sum(sin(pi*Foam::pow(f, m)));
    List<scalar> sol = A.solve();
    K1_ = sol[0];
    K2_ = sol[1];
}


void Foam::fv::LeishmanBeddoes::calcSeparated()
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
    if (tau_ > 0 and tau_ <= Tvl_)
    {
        Tf = 0.5*Tf_;
    }
    else if (tau_ > Tvl_ and tau_ <= 2.0*Tvl_)
    {
        Tf = 4.0*Tf_;
    }
    if (mag(alpha_) < mag(alphaPrev_) and mag(CNPrime_) < CN1_)
    {
        Tf = 0.5*Tf_;
    }

    // Calculate dynamic separation point
    DF_ = DFPrev_*exp(-deltaS_/Tf)
        + (fPrime_ - fPrimePrev_)*exp(-deltaS_/(2.0*Tf));
    fDoublePrime_ = fPrime_ - DF_;
    if (fDoublePrime_ < 0)
    {
        fDoublePrime_ = 0.0;
    }
    else if (fDoublePrime_ > 1)
    {
        fDoublePrime_ = 1.0;
    }

    // Calculate normal force coefficient including dynamic separation point
    CNF_ = CNAlpha_*alphaEquiv_*pow(((1.0 + sqrt(fDoublePrime_))/2.0), 2.0)
         + CNI_;

    if (debug)
    {
        Info<< "    DF: " << DF_ << endl;
        Info<< "    CNF: " << CNF_ << endl;
    }

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

    // Calculate Strouhal number time constant and set tau to zero to
    // allow multiple vortex shedding
    scalar Tst = 2.0*(1.0 - fDoublePrime_)/0.19;
    if (tau_ > (Tvl_ + Tst))
    {
        tau_ = 0.0;
    }

    // Evaluate vortex lift contributions, which are only increasing if angle
    // of attack increased in magnitude beyond a threshold
    scalar Tv = Tv_;
    if (tau_ < Tvl_ and (mag(alpha_) > mag(alphaPrev_)))
    {
        // Halve Tv if dAlpha/dt changes sign
        if (sign(deltaAlpha_) != sign(deltaAlphaPrev_))
        {
            Tv = 0.5*Tv_;
        }
        scalar KN = magSqr((1.0 + sqrt(fDoublePrime_)))/4.0;
        CV_ = CNC_*(1.0 - KN);
        CNV_ = CNVPrev_*exp(-deltaS_/Tv)
             + (CV_ - CVPrev_)*exp(-deltaS_/(2.0*Tv));
    }
    else
    {
        Tv = 0.5*Tv_;
        CV_ = 0.0;
        CNV_ = CNVPrev_*exp(-deltaS_/Tv);
    }

    // Total normal force coefficient is the combination of that from
    // circulatory effects, impulsive effects, dynamic separation, and vortex
    // lift
    CN_ = CNF_ + CNV_;

    // Calculate moment coefficient
    scalar pi = Foam::constant::mathematical::pi;
    scalar m = cmFitExponent_;
    scalar cmf = (K0_ + K1_*(1 - fDoublePrime_)
               + K2_*sin(pi*Foam::pow(fDoublePrime_, m)))*CNC_
               + profileData_.zeroLiftMomentCoeff();
    scalar cpv = 0.20*(1 - cos(pi*tau_/Tvl_));
    scalar cmv = -cpv*CNV_;
    CM_ = cmf + cmv;

    if (debug)
    {
        Info<< "    CV: " << CV_ << endl;
        Info<< "    CNV: " << CNV_ << endl;
    }
}


void Foam::fv::LeishmanBeddoes::update()
{
    timePrev_ = time_.value();
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
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::LeishmanBeddoes::LeishmanBeddoes
(
    const dictionary& dict,
    const word& modelName,
    const Time& time,
    profileData& profileData
)
:
    dynamicStallModel(dict, modelName, time, profileData),
    X_(0.0),
    XPrev_(0.0),
    Y_(0.0),
    YPrev_(0.0),
    A1_(coeffs_.lookupOrDefault("A1", 0.3)),
    A2_(coeffs_.lookupOrDefault("A2", 0.7)),
    b1_(coeffs_.lookupOrDefault("b1", 0.14)),
    b2_(coeffs_.lookupOrDefault("b2", 0.53)),
    alpha_(0.0),
    deltaAlpha_(0.0),
    deltaAlphaPrev_(0.0),
    a_(coeffs_.lookupOrDefault("speedOfSound", 1e12)),
    timePrev_(startTime_),
    D_(0.0),
    DPrev_(0.0),
    DP_(0.0),
    DPPrev_(0.0),
    CNP_(0.0),
    CNPPrev_(0.0),
    fPrime_(1.0),
    fPrimePrev_(1.0),
    DF_(0.0),
    DFPrev_(0.0),
    CV_(0.0),
    CVPrev_(0.0),
    CNV_(0.0),
    CNVPrev_(0.0),
    eta_(coeffs_.lookupOrDefault("eta", 0.95)),
    S1_(VGREAT),
    S2_(VGREAT),
    CD0_(VGREAT),
    stalledPrev_(false),
    Tp_(coeffs_.lookupOrDefault("Tp", 1.7)),
    Tf_(coeffs_.lookupOrDefault("Tf", 3.0)),
    Tv_(coeffs_.lookupOrDefault("Tv", 6.0)),
    Tvl_(coeffs_.lookupOrDefault("Tvl", 7.0)),
    tau_(0.0),
    tauPrev_(0.0),
    nNewTimes_(0),
    fCrit_(0.7),
    K0_(1e-6),
    K1_(0.0),
    K2_(0.0),
    cmFitExponent_(coeffs_.lookupOrDefault("cmFitExponent", 2)),
    CM_(0.0),
    Re_(0.0)
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
    scalar magU,
    scalar alphaDeg,
    scalar& cl,
    scalar& cd,
    scalar& cm
)
{
    scalar pi = Foam::constant::mathematical::pi;
    scalar time = time_.value();
    deltaT_ = time_.deltaT().value();

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
        alpha_ = alphaDeg/180.0*pi;
        alphaPrev_ = alpha_;
    }

    magU_ = magU;
    alpha_ = alphaDeg/180.0*pi;
    M_ = magU/a_;
    deltaAlpha_ = alpha_ - alphaPrev_;
    deltaS_ = 2*magU*deltaT_/c_;

    if (debug)
    {
        scalar cn0 = cl*cos(alpha_) - cd*sin(alpha_);
        Info<< "    -------------------------------------------------" << endl;
        Info<< "    " << modelName_ << " dynamic stall model correcting"
            << endl;
        Info<< "    -------------------------------------------------" << endl;
        Info<< "    New times: " << nNewTimes_ << endl;
        Info<< "    Time: " << time << endl;
        Info<< "    deltaT: " << deltaT_ << endl;
        Info<< "    deltaS: " << deltaS_ << endl;
        Info<< "    magU: " << magU << endl;
        Info<< "    Angle of attack (deg): " << alphaDeg << endl;
        Info<< "    alpha (rad): " << alpha_ << endl;
        Info<< "    deltaAlpha: " << deltaAlpha_ << endl;
        Info<< "    Mach number: " << M_ << endl;
        Info<< "    Initial normal force coefficient: " << cn0 << endl;
        Info<< "    Initial lift coefficient: " << cl << endl;
        Info<< "    Initial drag coefficient: " << cd << endl;
        Info<< "    Initial moment coefficient: " << cm << endl;
    }

    bool doCalcAlphaEquiv = coeffs_.lookupOrDefault("calcAlphaEquiv", false);
    if (doCalcAlphaEquiv)
    {
        calcAlphaEquiv();
    }
    else
    {
        alphaEquiv_ = alpha_;
    }
    // Evaluate static coefficient data if it has changed, e.g., from a
    // Reynolds number correction
    if (profileData_.staticStallAngleRad() != alphaSS_)
    {
        evalStaticData();
    }
    calcUnsteady();
    calcSeparated();

    // Modify coefficients
    cl = CN_*cos(alpha_) + CT_*sin(alpha_);
    cd = CN_*sin(alpha_) - CT_*cos(alpha_) + CD0_;
    cm = CM_;

    if (debug)
    {
        scalar alphE = alphaEquiv_/pi*180.0;
        Info<< "    Stalled: " << stalled_ << endl;
        Info<< "    tau: " << tau_ << endl;
        Info<< "    Equivalent angle of attack: " << alphE << endl;
        Info<< "    alphaPrime: " << alphaPrime_ << endl;
        Info<< "    fPrime: " << fPrime_ << endl;
        Info<< "    fDoublePrime: " << fDoublePrime_ << endl;
        Info<< "    Corrected normal force coefficient: " << CN_ << endl;
        Info<< "    Circulatory normal force coefficient: " << CNC_ << endl;
        Info<< "    Separation normal force coefficient: " << CNF_ << endl;
        Info<< "    Impulsive normal force coefficient: " << CNI_ << endl;
        Info<< "    Vortex normal force coefficient: " << CNV_ << endl;
        Info<< "    Tangential force coefficient: " << CT_ << endl;
        Info<< "    Corrected lift coefficient: " << cl << endl;
        Info<< "    Corrected drag coefficient: " << cd << endl;
        Info<< "    Corrected moment coefficient: " << cm << endl;
        Info<< "    -------------------------------------------------" << endl;
    }
}


void Foam::fv::LeishmanBeddoes::reduceParallel(bool inMesh)
{
    if (not inMesh)
    {
        alpha_ = VGREAT;
        alphaPrev_ = VGREAT;
        timePrev_ = VGREAT;
        X_ = VGREAT;
        XPrev_ = VGREAT;
        Y_ = VGREAT;
        YPrev_ = VGREAT;
        D_ = VGREAT;
        DPrev_ = VGREAT;
        DP_ = VGREAT;
        DPPrev_ = VGREAT;
        CNP_ = VGREAT;
        CNPPrev_ = VGREAT;
        DF_ = VGREAT;
        DFPrev_ = VGREAT;
        fPrime_ = VGREAT;
        fPrimePrev_ = VGREAT;
        CV_ = VGREAT;
        CVPrev_ = VGREAT;
        CNV_ = VGREAT;
        CNVPrev_ = VGREAT;
        stalled_ = false;
        stalledPrev_ = false;
        tau_ = VGREAT;
        tauPrev_ = VGREAT;
    }

    reduce(alpha_, minOp<scalar>());
    reduce(alphaPrev_, minOp<scalar>());
    reduce(timePrev_, minOp<scalar>());
    reduce(X_, minOp<scalar>());
    reduce(XPrev_, minOp<scalar>());
    reduce(Y_, minOp<scalar>());
    reduce(YPrev_, minOp<scalar>());
    reduce(D_, minOp<scalar>());
    reduce(DPrev_, minOp<scalar>());
    reduce(DP_, minOp<scalar>());
    reduce(DPPrev_, minOp<scalar>());
    reduce(CNP_, minOp<scalar>());
    reduce(CNPPrev_, minOp<scalar>());
    reduce(DF_, minOp<scalar>());
    reduce(DFPrev_, minOp<scalar>());
    reduce(fPrime_, minOp<scalar>());
    reduce(fPrimePrev_, minOp<scalar>());
    reduce(CV_, minOp<scalar>());
    reduce(CVPrev_, minOp<scalar>());
    reduce(CNV_, minOp<scalar>());
    reduce(CNVPrev_, minOp<scalar>());
    reduce(stalled_, maxOp<bool>());
    reduce(stalledPrev_, maxOp<bool>());
    reduce(tau_, minOp<scalar>());
    reduce(tauPrev_, minOp<scalar>());
}

// ************************************************************************* //
