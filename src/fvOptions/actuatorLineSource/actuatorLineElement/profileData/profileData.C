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

#include "profileData.H"
#include "interpolateUtils.H"
#include "simpleMatrix.H"
#include "unitConversion.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::profileData::interpolate
(
    scalar xNew,
    List<scalar>& xOld,
    List<scalar>& yOld
)
{
    label index = 0;
    label indexP = 0;
    label indexM = 0;
    scalar error = 1.0E30;
    forAll(xOld, i)
    {
        scalar diff = mag(xNew - xOld[i]);
        if (diff < error)
        {
            index = i;
            error = diff;
        }
    }
    if (xNew < xOld[index])
    {
        if (index == 0)
        {
            indexP = 1;
            indexM = indexP - 1;
        }
        else
        {
            indexP = index;
            indexM = indexP - 1;
        }
        return yOld[indexM]
               + ((yOld[indexP]
               - yOld[indexM])/(xOld[indexP]
               - xOld[indexM]))*(xNew - xOld[indexM]);
    }
    else if (xNew > xOld[index])
    {
        if (index == xOld.size() - 1)
        {
            indexP = xOld.size() - 1;
            indexM = indexP - 1;
        }
        else
        {
            indexP = index + 1;
            indexM = indexP - 1;
        }
        return yOld[indexM] + ((yOld[indexP]
               - yOld[indexM])/(xOld[indexP]
               - xOld[indexM]))*(xNew - xOld[indexM]);
    }
    else if (xNew == xOld[index])
    {
        return yOld[index];
    }
    else
    {
        return 0.0;
    }
}


void Foam::profileData::read2DArray
(
    List<List<scalar> >& data,
    word keyword
)
{
    List<List<scalar> > inputData = dict_.lookup(keyword);
    // Check that the number of columns is correct based on the Re list
    if ((inputData[0].size() - 1) != ReList_.size())
    {
        FatalErrorIn("void profileData::read2DArray()")
            << "Number of columns in " << keyword
            << " does not match number of Reynolds numbers"
            << abort(FatalError);
    }

    data.setSize(inputData.size());
    forAll(data, i)
    {
        // Data should not include first column
        data[i].setSize(inputData[i].size() - 1);
        forAll(data[i], j)
        {
            data[i][j] = inputData[i][j+1];
        }
    }
}


void Foam::profileData::readSingleRe()
{
    // Read reference Reynolds number, and if present turn on Reynolds number
    // corrections
    ReRef_ = dict_.lookupOrDefault("Re", VSMALL);
    Re_ = ReRef_;
    correctRe_ = (ReRef_ > VSMALL);

    // Look up sectional coefficient data
    List<List<scalar> > coefficientData = dict_.lookup("data");

    // Create lists from coefficient data
    angleOfAttackListOrg_.setSize(coefficientData.size());
    liftCoefficientListOrg_.setSize(coefficientData.size());
    dragCoefficientListOrg_.setSize(coefficientData.size());
    momentCoefficientListOrg_.setSize(coefficientData.size());
    forAll(coefficientData, i)
    {
        angleOfAttackListOrg_[i] = coefficientData[i][0];
        liftCoefficientListOrg_[i] = coefficientData[i][1];
        dragCoefficientListOrg_[i] = coefficientData[i][2];

        if (coefficientData[i].size() > 3)
        {
            momentCoefficientListOrg_[i] = coefficientData[i][3];
        }
        else
        {
            momentCoefficientListOrg_[i] = VSMALL;
        }
    }

    // Initially lists are identical to original
    angleOfAttackList_ = angleOfAttackListOrg_;
    liftCoefficientList_ = liftCoefficientListOrg_;
    dragCoefficientList_ = dragCoefficientListOrg_;
    momentCoefficientList_ = momentCoefficientListOrg_;
}


void Foam::profileData::readMultiRe()
{
    // Turn off Reynolds number corrections since these will be interpolated
    correctRe_ = false;

    // Read list of Reynolds numbers for dataset
    dict_.lookup("ReList") >> ReList_;

    // Scale Reynolds numbers if desired
    // This functionality can be used to approximate full scale loading without
    // changing the flow Reynolds number
    // Lookup table will effectively use Reynolds numbers ReScale times
    // larger than the computed element Reynolds number
    if (dict_.found("ReScale"))
    {
        scalar ReScale = 1.0;
        dict_.lookup("ReScale") >> ReScale;
        forAll(ReList_, i)
        {
            ReList_[i] = ReList_[i] / ReScale;
        }
    }

    // Define keywords for coefficient data
    word clKeyword = "clData";
    word cdKeyword = "cdData";
    word cmKeyword = "cmData";

    // Create angle of attack lists, and ensure these are equal for all
    // quantities
    angleOfAttackList_ = readAngleOfAttackList(clKeyword);
    if (angleOfAttackList_ != readAngleOfAttackList(cdKeyword))
    {
        FatalErrorIn("void profileData::readMultiRe()")
            << "Angle of attack lists do not match"
            << abort(FatalError);
    }

    // Set sizes of coefficient lists
    liftCoefficientList_.setSize(angleOfAttackList_.size());
    dragCoefficientList_.setSize(angleOfAttackList_.size());
    momentCoefficientList_.setSize(angleOfAttackList_.size());

    // Look up 2-D array of lift coefficient data
    read2DArray
    (
        liftCoefficientLists_,
        clKeyword
    );

    // Look up 2-D array of drag coefficient data
    read2DArray
    (
        dragCoefficientLists_,
        cdKeyword
    );

    if (liftCoefficientLists_.size() != dragCoefficientLists_.size())
    {
        FatalErrorIn("void profileData::read()")
            << "Lift and drag coefficient data must be the same size"
            << abort(FatalError);
    }

    // Look up 2-D array of moment coefficient data if available
    if (dict_.found(cmKeyword))
    {
        read2DArray
        (
            momentCoefficientLists_,
            cmKeyword
        );
    }
    else
    {
        momentCoefficientLists_ = liftCoefficientLists_;
        forAll(momentCoefficientLists_, i)
        {
            momentCoefficientLists_[i] = liftCoefficientLists_[i]*0.0;
        }
    }

    // Calculate static stall angle, zero lift AoA, etc. for all Re
    analyzeMultiRe();
}


Foam::List<scalar> Foam::profileData::readAngleOfAttackList(word keyword)
{
    List<List<scalar> > data = dict_.lookup(keyword);
    List<scalar> alphaList(data.size());
    forAll(alphaList, i)
    {
        alphaList[i] = data[i][0];
    }
    return alphaList;
}


void Foam::profileData::calcStaticStallAngle()
{
    // Static stall is where the slope of the drag coefficient curve first
    // breaks {threshold} per degree
    scalar threshold = 0.03;
    scalar alpha=GREAT, cd0, cd1, slope, dAlpha;
    forAll(angleOfAttackList_, i)
    {
        alpha = angleOfAttackList_[i];
        if (alpha > 2 && alpha < 30)
        {
            cd1 = dragCoefficient(alpha + 1.0);
            cd0 = dragCoefficient(alpha);
            dAlpha = 1.0;
            slope = (cd1 - cd0)/dAlpha;
            if (slope > threshold)
            {
                staticStallAngle_ = alpha;
                break;
            }
        }
    }
}


void Foam::profileData::calcZeroLiftDragCoeff()
{
    List<scalar> cdList = dragCoefficientList(-10, 10);
    List<scalar> clList = liftCoefficientList(-10, 10);
    zeroLiftDragCoeff_ = interpolate(0, clList, cdList);
}


void Foam::profileData::calcZeroLiftAngleOfAttack()
{
    List<scalar> clList = liftCoefficientList(-10, 10);
    List<scalar> alphaList = angleOfAttackList(-10, 10);
    zeroLiftAngleOfAttack_ = interpolate(0, clList, alphaList);
}


void Foam::profileData::calcZeroLiftMomentCoeff()
{
    List<scalar> cmList = momentCoefficientList(-10, 10);
    List<scalar> clList = liftCoefficientList(-10, 10);
    zeroLiftMomentCoeff_ = interpolate(0, clList, cmList);
}


void Foam::profileData::calcNormalCoeffSlope()
{
    scalar alphaHigh = 0.5*staticStallAngle_;
    List<scalar> alphaList = degToRad(angleOfAttackList(0, alphaHigh));
    List<scalar> cnList = normalCoefficientList(0, alphaHigh);
    simpleMatrix<scalar> A(2);
    A[0][0] = alphaList.size();
    A[0][1] = Foam::sum(alphaList);
    A[1][0] = Foam::sum(alphaList);
    A[1][1] = Foam::sum(Foam::sqr(alphaList));
    A.source()[0] = Foam::sum(cnList);
    A.source()[1] = Foam::sum(cnList*alphaList);
    normalCoeffSlope_ = A.solve()[1];
}


void Foam::profileData::interpCoeffLists()
{
    // Create lists for lift, drag, and moment coefficient at current Re
    forAll(angleOfAttackList_, i)
    {
        liftCoefficientList_[i] = interpolateUtils::interpolate1D
        (
            Re_,
            ReList_,
            liftCoefficientLists_[i]
        );
        dragCoefficientList_[i] = interpolateUtils::interpolate1D
        (
            Re_,
            ReList_,
            dragCoefficientLists_[i]
        );
        momentCoefficientList_[i] = interpolateUtils::interpolate1D
        (
            Re_,
            ReList_,
            momentCoefficientLists_[i]
        );
    }
}


void Foam::profileData::interpPropsMultiRe()
{
    // Since Reynolds number list is the same for all interpolations,
    // precompute interpolation data for faster interpolation
    label interpIndex = interpolateUtils::binarySearch(ReList_, Re_);
    scalar interpFraction = interpolateUtils::getPart
    (
        Re_,
        ReList_,
        interpIndex
    );

    staticStallAngle_ = interpolateUtils::interpolate1D
    (
        interpFraction,
        staticStallAngleList_,
        interpIndex
    );
    zeroLiftDragCoeff_ = interpolateUtils::interpolate1D
    (
        interpFraction,
        zeroLiftDragCoeffList_,
        interpIndex
    );
    zeroLiftAngleOfAttack_ = interpolateUtils::interpolate1D
    (
        interpFraction,
        zeroLiftAngleOfAttackList_,
        interpIndex
    );
    zeroLiftMomentCoeff_ = interpolateUtils::interpolate1D
    (
        interpFraction,
        zeroLiftMomentCoeffList_,
        interpIndex
    );
    normalCoeffSlope_ = interpolateUtils::interpolate1D
    (
        interpFraction,
        normalCoeffSlopeList_,
        interpIndex
    );
}


void Foam::profileData::analyzeMultiRe()
{
    scalar ReOld = Re_;
    forAll(ReList_, i)
    {
        Re_ = ReList_[i];
        interpCoeffLists();
        analyze();
        staticStallAngleList_.append(staticStallAngle_);
        zeroLiftDragCoeffList_.append(zeroLiftDragCoeff_);
        zeroLiftAngleOfAttackList_.append(zeroLiftAngleOfAttack_);
        zeroLiftMomentCoeffList_.append(zeroLiftMomentCoeff_);
        normalCoeffSlopeList_.append(normalCoeffSlope_);
    }
    Re_ = ReOld;
}


Foam::List<scalar> Foam::profileData::subList
(
    scalar alphaDegStart,
    scalar alphaDegStop,
    const List<scalar>& fullList
)
{
    List<scalar> newList;
    forAll(angleOfAttackList_, i)
    {
        if
        (
            angleOfAttackList_[i] >= alphaDegStart
            and
            angleOfAttackList_[i] <= alphaDegStop
        )
        {
            newList.append(fullList[i]);
        }
    }
    return newList;
}


Foam::scalar Foam::profileData::convertToCN
(
    scalar cl,
    scalar cd,
    scalar angleOfAttackDeg
)
{
    scalar angleOfAttackRad = degToRad(angleOfAttackDeg);
    return cl*cos(angleOfAttackRad) + cd*sin(angleOfAttackRad);
}


Foam::scalar Foam::profileData::convertToCC
(
    scalar cl,
    scalar cd,
    scalar angleOfAttackDeg
)
{
    scalar angleOfAttackRad = degToRad(angleOfAttackDeg);
    return cl*sin(angleOfAttackRad) - cd*cos(angleOfAttackRad);
}


Foam::scalar Foam::profileData::convertToCL
(
    scalar cn,
    scalar cc,
    scalar angleOfAttackDeg
)
{
    scalar angleOfAttackRad = degToRad(angleOfAttackDeg);
    return cn*cos(angleOfAttackRad) + cc*sin(angleOfAttackRad);
}


Foam::scalar Foam::profileData::convertToCD
(
    scalar cn,
    scalar cc,
    scalar angleOfAttackDeg
)
{
    scalar angleOfAttackRad = degToRad(angleOfAttackDeg);
    return cn*sin(angleOfAttackRad) - cc*cos(angleOfAttackRad);
}


Foam::scalar Foam::profileData::convertToCRT
(
    scalar cl,
    scalar cd,
    scalar inflowVelAngleDeg
)
{
    scalar inflowVelAngleRad = degToRad(inflowVelAngleDeg);
    return cl*sin(inflowVelAngleRad) - cd*cos(inflowVelAngleRad);
}


Foam::scalar Foam::profileData::convertToCRN
(
    scalar cl,
    scalar cd,
    scalar inflowVelAngleDeg
)
{
    scalar inflowVelAngleRad = degToRad(inflowVelAngleDeg);
    return cl*cos(inflowVelAngleRad) + cd*sin(inflowVelAngleRad);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::profileData::profileData
(
    const word& name,
    const dictionary& dict,
    const label& debug
)
:
    name_(name),
    dict_(dict),
    debug(debug),
    tableType_(dict.lookupOrDefault("tableType", word("singleRe"))),
    Re_(VSMALL),
    ReRef_(VSMALL),
    correctRe_(false),
    staticStallAngle_(VGREAT),
    zeroLiftDragCoeff_(VGREAT),
    zeroLiftAngleOfAttack_(VGREAT),
    zeroLiftMomentCoeff_(VGREAT),
    normalCoeffSlope_(VGREAT)
{
    if (tableType_ == "singleRe")
    {
        readSingleRe();
    }
    else if (tableType_ == "multiRe")
    {
        readMultiRe();
    }
    else
    {
        FatalErrorIn("Foam::profileData::profileData")
            << "Unknown profileData tableType " << tableType_
            << ". Must be either 'singleRe' or 'multiRe'"
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::profileData>
Foam::profileData::New
(
    const word& name,
    const dictionary& dict,
    const label& debug
)
{
    return autoPtr<profileData>
    (
        new profileData
        (
            name,
            dict,
            debug
        )
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::profileData::~profileData()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::profileData::analyze()
{
    calcStaticStallAngle();
    calcZeroLiftDragCoeff();
    calcZeroLiftAngleOfAttack();
    calcZeroLiftMomentCoeff();
    calcNormalCoeffSlope();
}


Foam::scalar Foam::profileData::liftCoefficient(scalar angleOfAttackDeg)
{
    return interpolate
    (
        angleOfAttackDeg,
        angleOfAttackList_,
        liftCoefficientList_
    );
}


Foam::scalar Foam::profileData::dragCoefficient(scalar angleOfAttackDeg)
{
    return interpolate
    (
        angleOfAttackDeg,
        angleOfAttackList_,
        dragCoefficientList_
    );
}


Foam::scalar Foam::profileData::momentCoefficient(scalar angleOfAttackDeg)
{
    return interpolate
    (
        angleOfAttackDeg,
        angleOfAttackList_,
        momentCoefficientList_
    );
}


Foam::scalar Foam::profileData::normalCoefficient(scalar angleOfAttackDeg)
{
    return convertToCN
    (
       liftCoefficient(angleOfAttackDeg),
       dragCoefficient(angleOfAttackDeg),
       angleOfAttackDeg
    );
}


Foam::scalar Foam::profileData::chordwiseCoefficient(scalar angleOfAttackDeg)
{
    return convertToCC
    (
       liftCoefficient(angleOfAttackDeg),
       dragCoefficient(angleOfAttackDeg),
       angleOfAttackDeg
    );
}


void Foam::profileData::updateRe(scalar Re)
{
    if (tableType_ == "singleRe" and correctRe_ and Re != Re_)
    {
        Re_ = Re;

        // Correct drag coefficients
        scalar fReRef = Foam::pow((Foam::log(ReRef_) - 0.407), -2.64);
        scalar fRe = Foam::pow((Foam::log(Re) - 0.407), -2.64);
        scalar K = fReRef/fRe;
        dragCoefficientList_ = dragCoefficientListOrg_/K;

        if (debug)
        {
            Info<< "Correcting " << name_ << " profile data for Reynolds number"
                << " effects" << endl;
            Info<< "    Re: " << Re << endl;
            Info<< "    ReRef: " << ReRef_ << endl;
            Info<< "    f(Re): " << fRe << endl;
            Info<< "    f(ReRef): " << fReRef << endl;
            Info<< "    K (drag): " << K << endl;
        }

        // Correct lift coefficients
        scalar n = dict_.lookupOrDefault("liftReCorrExp", 0.23);
        K = pow((Re/ReRef_), n);
        forAll(liftCoefficientList_, i)
        {
            scalar alphaNew = angleOfAttackListOrg_[i]/K;
            liftCoefficientList_[i] = interpolate
            (
                alphaNew,
                angleOfAttackListOrg_,
                liftCoefficientListOrg_
            );
            liftCoefficientList_[i] *= K;
        }

        if (debug)
        {
            Info<< "    n: " << n << endl;
            Info<< "    K (lift): " << K << endl;
            Info<< "    Initial minimum drag coefficient: "
                << Foam::min(dragCoefficientListOrg_) << endl;
            Info<< "    Corrected minimum drag coefficient: "
                << Foam::min(dragCoefficientList_) << endl;
            Info<< "    Initial maximum lift coefficient: "
                << Foam::max(liftCoefficientListOrg_) << endl;
            Info<< "    Corrected maximum lift coefficient: "
                << Foam::max(liftCoefficientList_) << endl;
        }

        // Recalculate static stall angle, etc.
        analyze();
    }
    else if (tableType_ == "multiRe" and Re != Re_)
    {
        Re_ = Re;
        interpCoeffLists();
        interpPropsMultiRe();
    }
    else
    {
        Re_ = Re;
    }
}


const Foam::dictionary& Foam::profileData::dict()
{
    return dict_;
}


const Foam::List<scalar>& Foam::profileData::angleOfAttackList()
{
    return angleOfAttackList_;
}


const Foam::List<scalar>& Foam::profileData::liftCoefficientList()
{
    return liftCoefficientList_;
}


const Foam::List<scalar>& Foam::profileData::dragCoefficientList()
{
    return dragCoefficientList_;
}


const Foam::List<scalar>& Foam::profileData::momentCoefficientList()
{
    return momentCoefficientList_;
}


Foam::List<scalar> Foam::profileData::angleOfAttackList
(
    scalar alphaDegStart,
    scalar alphaDegStop
)
{
    return subList
    (
        alphaDegStart,
        alphaDegStop,
        angleOfAttackList_
    );
}


Foam::List<scalar> Foam::profileData::liftCoefficientList
(
    scalar alphaDegStart,
    scalar alphaDegStop
)
{
    return subList
    (
        alphaDegStart,
        alphaDegStop,
        liftCoefficientList_
    );
}


Foam::List<scalar> Foam::profileData::dragCoefficientList
(
    scalar alphaDegStart,
    scalar alphaDegStop
)
{
    return subList
    (
        alphaDegStart,
        alphaDegStop,
        dragCoefficientList_
    );
}


Foam::List<scalar> Foam::profileData::momentCoefficientList
(
    scalar alphaDegStart,
    scalar alphaDegStop
)
{
    return subList
    (
        alphaDegStart,
        alphaDegStop,
        momentCoefficientList_
    );
}


Foam::List<scalar> Foam::profileData::normalCoefficientList
(
    scalar alphaDegStart,
    scalar alphaDegStop
)
{
    List<scalar> cnList(liftCoefficientList_.size());
    forAll(cnList, i)
    {
        cnList[i] = convertToCN
        (
            liftCoefficientList_[i],
            dragCoefficientList_[i],
            angleOfAttackList_[i]
        );
    }
    return subList
    (
        alphaDegStart,
        alphaDegStop,
        cnList
    );
}


Foam::List<scalar> Foam::profileData::chordwiseCoefficientList
(
    scalar alphaDegStart,
    scalar alphaDegStop
)
{
    List<scalar> ccList(liftCoefficientList_.size());
    forAll(ccList, i)
    {
        ccList[i] = convertToCC
        (
            liftCoefficientList_[i],
            dragCoefficientList_[i],
            angleOfAttackList_[i]
        );
    }
    return subList
    (
        alphaDegStart,
        alphaDegStop,
        ccList
    );
}


Foam::scalar Foam::profileData::staticStallAngleRad()
{
    if (staticStallAngle_ == VGREAT)
    {
        calcStaticStallAngle();
    }
    return degToRad(staticStallAngle_);
}


Foam::scalar Foam::profileData::zeroLiftDragCoeff()
{
    if (zeroLiftDragCoeff_ == VGREAT)
    {
        calcZeroLiftDragCoeff();
    }
    return zeroLiftDragCoeff_;
}


Foam::scalar Foam::profileData::zeroLiftAngleOfAttack()
{
    if (zeroLiftAngleOfAttack_ == VGREAT)
    {
        calcZeroLiftAngleOfAttack();
    }
    return zeroLiftAngleOfAttack_;
}


Foam::scalar Foam::profileData::zeroLiftMomentCoeff()
{
    if (zeroLiftMomentCoeff_ == VGREAT)
    {
        calcZeroLiftMomentCoeff();
    }
    return zeroLiftMomentCoeff_;
}


Foam::scalar Foam::profileData::normalCoeffSlope()
{
    if (normalCoeffSlope_ == VGREAT)
    {
        calcNormalCoeffSlope();
    }
    return normalCoeffSlope_;
}


Foam::scalar Foam::profileData::Re()
{
    return Re_;
}


bool Foam::profileData::correctRe()
{
    return correctRe_;
}


// ************************************************************************* //
