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


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

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
        if(diff < error)
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


/*reads matrix data from file, data is assumed to have index values in 
first row and column, where the row will be read to xvalues and the column
to y values, and data is given as data[y][x]*/
void Foam::profileData::readMatrix
(
    List<scalar> &xvalues,
    List<scalar> &yvalues,
    List<List<scalar> > &data,
    word keyword
)
{
    bool buildXList = xvalues.size() == 0;
    bool buildYList = yvalues.size() == 0;
    List<List<scalar> > empty;
    List<List<scalar> > coefficientData = 
        dict_.lookupOrDefault<List<List<scalar> > >(keyword, empty);

    if (coefficientData.size() > 0)
    {
        //for simplicity, we use the same index lists for all coefficients, 
        //hence if the lists are constructed
        if (buildXList)
            yvalues.setSize(coefficientData.size() - 1);
        if (buildYList)
            xvalues.setSize(coefficientData[0].size() - 1);
        data.setSize(coefficientData.size() - 1);
        for (int i = 1; i < coefficientData[0].size(); i++)
        {
            if(buildXList)
                xvalues[i-1] = coefficientData[0][i];
            else
            {
                if(xvalues[i-1] != coefficientData[0][i])
                {
                    word errorMessage = 
                        string("Index elements in ") + keyword + 
                        " must be identical to other coefficient lists";
                    error myerror(errorMessage);
                    myerror.abort();
                }
            }
            if (i > 1 && xvalues[i-1] < xvalues[i-2])
            {
                std::string errorMessage = 
                    string("Index elements in ") + keyword + 
                    " must be ordered with smallest element first";
                error myerror(errorMessage);
                myerror.abort();
            }
        }
        for (int i = 1; i < coefficientData.size(); i++)
        {
            if (buildYList)
                yvalues[i-1] = coefficientData[i][0];
            else
            {
                if(yvalues[i-1] != coefficientData[i][0])
                {
                    std::string errorMessage = 
                          string("Index elements in ") + keyword + 
                          " must be identical to other coefficient lists";
                    error myerror(errorMessage);
                    myerror.abort();
                }
            }
            if (i > 1 && yvalues[i-1] < yvalues[i-2])
            {
                std::string errorMessage = 
                    string("Index elements in ") + keyword + 
                    " must be ordered with smallest element first";
                error myerror(errorMessage);
                myerror.abort();
            }
            data[i-1].setSize(coefficientData[i-1].size() - 1);
            if (coefficientData[i-1].size() != coefficientData[0].size())
            {
                std::string errorMessage = 
                    string("Element size of data in ") + keyword + 
                    " varies in size, all elements must have equal size";
                error myerror(errorMessage);
                myerror.abort();
            }
            for (int j = 1; j < coefficientData[i].size(); j++)
            {
                data[i-1][j-1] = coefficientData[i][j];
            }
        }
    }
}

void Foam::profileData::read()
{
    // Read reference Reynolds number, and if present turn on Reynolds number
    // corrections
    ReRef_ = dict_.lookupOrDefault("Re", VSMALL);
    correctRe_ = (ReRef_ > VSMALL);
    
    //Look up matrix data for Cl
    readMatrix
    (
        ReynoldsNumberListMatrixOrg_, 
        angleOfAttackListMatrixOrg_, 
        liftCoefficientMatrixOrg_, 
        "dataCl"
    );

    //Look up matrix data for Cd
    readMatrix
    (
        ReynoldsNumberListMatrixOrg_, 
        angleOfAttackListMatrixOrg_, 
        dragCoefficientMatrixOrg_, 
        "dataCd"
    );

    if (liftCoefficientMatrixOrg_.size() != dragCoefficientMatrixOrg_.size())
    {
        error myerror
        (
            "if lift and drag data is given in matrix form, "
            "both lift and drag data has to be given in this form"
        );
        myerror.abort();
    }
    
    //Look up matrix data for Moment
    readMatrix
    (
        ReynoldsNumberListMatrixOrg_, 
        angleOfAttackListMatrixOrg_, 
        momentCoefficientMatrixOrg_, 
        "dataMoment"
    );
    
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

    if (liftCoefficientMatrixOrg_.size())
    {
        buildReynoldsList();
    }
}


void Foam::profileData::calcStaticStallAngle()
{
    // Static stall is where the slope of the drag coefficient curve first
    // breaks {threshold} per degree
    scalar threshold = 0.03;
    scalar alpha=GREAT, cd0, cd1, slope, dAlpha;
    List<scalar>* angleOfAttackListptr;
    if (dragCoefficientMatrixOrg_.size() > 0)
        angleOfAttackListptr = &angleOfAttackListMatrixOrg_;
    else
        angleOfAttackListptr = &angleOfAttackList_;
    List<scalar>& angleOfAttackList = *angleOfAttackListptr;
    forAll(angleOfAttackList, i)
    {
        alpha = angleOfAttackList[i];
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


void Foam::profileData::calculateCoefficients()
{
    calcStaticStallAngle();
    calcZeroLiftDragCoeff();
    calcZeroLiftAngleOfAttack();
    calcZeroLiftMomentCoeff();
    calcNormalCoeffSlope();
}


void Foam::profileData::getInterpolatedCoefficients()
{
    staticStallAngle_ = interpolateUtils::interpolate1d
    (
        Re_, 
        ReynoldsNumberListMatrixOrg_, 
        staticStallAngleList_
    );
    zeroLiftDragCoeff_ = interpolateUtils::interpolate1d
    (
        Re_, 
        ReynoldsNumberListMatrixOrg_, 
        zeroLiftDragCoeffList_
    );
    zeroLiftAngleOfAttack_ = interpolateUtils::interpolate1d
    (
        Re_, 
        ReynoldsNumberListMatrixOrg_, 
        zeroLiftAngleOfAttackList_
    );
    zeroLiftMomentCoeff_ = interpolateUtils::interpolate1d
    (
        Re_, 
        ReynoldsNumberListMatrixOrg_, 
        zeroLiftMomentCoeffList_
    );
    normalCoeffSlope_ = interpolateUtils::interpolate1d
    (
        Re_, 
        ReynoldsNumberListMatrixOrg_, 
        normalCoeffSlopeList_
    );
}


void Foam::profileData::buildReynoldsList()
{
    scalar Re_old = Re_;
    forAll(ReynoldsNumberListMatrixOrg_, i)
    {
        Re_ = ReynoldsNumberListMatrixOrg_[i];
        calculateCoefficients();
        staticStallAngleList_.append(staticStallAngle_);
        zeroLiftDragCoeffList_.append(zeroLiftDragCoeff_);
        zeroLiftAngleOfAttackList_.append(zeroLiftAngleOfAttack_);
        zeroLiftMomentCoeffList_.append(zeroLiftMomentCoeff_);
        normalCoeffSlopeList_.append(normalCoeffSlope_);
    }
    Re_ = Re_old;
}


Foam::List<scalar> Foam::profileData::subList
(
    scalar alphaDegStart,
    scalar alphaDegStop,
    List<scalar>& oldList,
    List<List<scalar> >& oldMatrix
)
{
    List<scalar> newList;
    if (oldMatrix.size() > 0)
    {
        forAll(angleOfAttackListMatrixOrg_, i)
        {
            if 
            (
                angleOfAttackListMatrixOrg_[i] >= alphaDegStart
                and
                angleOfAttackListMatrixOrg_[i] <= alphaDegStop
            )
            {
                newList.append
                (
                    interpolateUtils::interpolate1d
                    (
                        Re_, 
                        ReynoldsNumberListMatrixOrg_, 
                        oldMatrix[i]
                    )
                );
            }
        }
    }
    else
    {
        forAll(angleOfAttackList_, i)
        {
            if 
            (
                angleOfAttackList_[i] >= alphaDegStart
                and
                angleOfAttackList_[i] <= alphaDegStop
            )
            {
                newList.append(oldList[i]);
            }
        }
    }
    return newList;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


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
    Re_(VSMALL),
    ReRef_(VSMALL),
    correctRe_(false),
    staticStallAngle_(VGREAT),
    zeroLiftDragCoeff_(VGREAT),
    zeroLiftAngleOfAttack_(VGREAT),
    zeroLiftMomentCoeff_(VGREAT),
    normalCoeffSlope_(VGREAT)
{
    read();
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
    if (staticStallAngleList_.size() > 0) //if interpolated tables exist
        getInterpolatedCoefficients();
    else
        calculateCoefficients();
}


Foam::scalar Foam::profileData::liftCoefficient(scalar angleOfAttackDeg)
{
    if (liftCoefficientMatrixOrg_.size() > 0)
    {
        return interpolateUtils::interpolate2d
        (
            Re_,
            angleOfAttackDeg,
            ReynoldsNumberListMatrixOrg_,
            angleOfAttackListMatrixOrg_,
            liftCoefficientMatrixOrg_
        );
    }
    else
    {
        return interpolate
        (
            angleOfAttackDeg,
            angleOfAttackList_,
            liftCoefficientList_
        );
    }
}


Foam::scalar Foam::profileData::dragCoefficient(scalar angleOfAttackDeg)
{
    if (dragCoefficientMatrixOrg_.size() > 0)
    {
        return interpolateUtils::interpolate2d
        (
            Re_,
            angleOfAttackDeg,
            ReynoldsNumberListMatrixOrg_,
            angleOfAttackListMatrixOrg_,
            dragCoefficientMatrixOrg_
        );
    }
    else
    {
        return interpolate
        (
            angleOfAttackDeg,
            angleOfAttackList_,
            dragCoefficientList_
        );
    }
}


Foam::scalar Foam::profileData::momentCoefficient(scalar angleOfAttackDeg)
{
    if (momentCoefficientMatrixOrg_.size() > 0)
    {
        return interpolateUtils::interpolate2d
        (
            Re_,
            angleOfAttackDeg,
            ReynoldsNumberListMatrixOrg_,
            angleOfAttackListMatrixOrg_,
            momentCoefficientMatrixOrg_
        );
    }
    else
    {
        return interpolate
        (
            angleOfAttackDeg,
            angleOfAttackList_,
            momentCoefficientList_
        );
    }
}


void Foam::profileData::updateRe(scalar Re)
{
    if (correctRe_ and Re != Re_)
    {
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

    Re_ = Re;
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
    List<scalar> newList;
    if (liftCoefficientMatrixOrg_.size() > 0)
    {
        forAll(angleOfAttackListMatrixOrg_, i)
        {
            if 
            (
                angleOfAttackListMatrixOrg_[i] >= alphaDegStart
                and
                angleOfAttackListMatrixOrg_[i] <= alphaDegStop
            )
            {
                newList.append(angleOfAttackListMatrixOrg_[i]);
            }
        }
    }
    else
    {
        forAll(angleOfAttackList_, i)
        {
            if 
            (
                angleOfAttackList_[i] >= alphaDegStart
                and
                angleOfAttackList_[i] <= alphaDegStop
            )
            {
                newList.append(angleOfAttackList_[i]);
            }
        }
    }
    return newList;
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
        liftCoefficientList_,
        liftCoefficientMatrixOrg_
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
        dragCoefficientList_,
        dragCoefficientMatrixOrg_
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
        momentCoefficientList_,
        momentCoefficientMatrixOrg_
    );
}


Foam::List<scalar> Foam::profileData::normalCoefficientList
(
    scalar alphaDegStart,
    scalar alphaDegStop
)
{
    List<scalar> newList;
    if (liftCoefficientMatrixOrg_.size() > 0)
    {
        forAll(angleOfAttackListMatrixOrg_, i)
        {
            if 
            (
                angleOfAttackListMatrixOrg_[i] >= alphaDegStart
                and
                angleOfAttackListMatrixOrg_[i] <= alphaDegStop
            )
            {
                newList.append
                (
                    interpolateUtils::interpolate1d
                    (
                        Re_, 
                        ReynoldsNumberListMatrixOrg_, 
                        liftCoefficientMatrixOrg_[i])*
                        cos(degToRad(angleOfAttackList_[i])
                    )
                  + interpolateUtils::interpolate1d
                    (
                        Re_, 
                        ReynoldsNumberListMatrixOrg_, 
                        dragCoefficientMatrixOrg_[i])*
                        sin(degToRad(angleOfAttackList_[i])
                    )
                );
            }
        }
    }
    else
    {
        forAll(angleOfAttackList_, i)
        {
            if 
            (
                angleOfAttackList_[i] >= alphaDegStart
                and
                angleOfAttackList_[i] <= alphaDegStop
            )
            {
                newList.append
                (
                    liftCoefficientList_[i]*cos(degToRad(angleOfAttackList_[i]))
                  + dragCoefficientList_[i]*sin(degToRad(angleOfAttackList_[i]))
                );
            }
        }
    }
    return newList;
}


Foam::List<scalar> Foam::profileData::chordwiseCoefficientList
(
    scalar alphaDegStart,
    scalar alphaDegStop
)
{
    List<scalar> newList;
    if(liftCoefficientMatrixOrg_.size() > 0)
    {
        forAll(angleOfAttackList_, i)
        {
            if 
            (
                angleOfAttackListMatrixOrg_[i] >= alphaDegStart
                and
                angleOfAttackListMatrixOrg_[i] <= alphaDegStop
            )
            {
                newList.append
                (
                    interpolateUtils::interpolate1d
                    (
                        Re_, 
                        ReynoldsNumberListMatrixOrg_, 
                        liftCoefficientMatrixOrg_[i])*
                        sin(degToRad(angleOfAttackList_[i])
                    )
                  - interpolateUtils::interpolate1d
                    (
                        Re_, 
                        ReynoldsNumberListMatrixOrg_, 
                        dragCoefficientMatrixOrg_[i])*
                        cos(degToRad(angleOfAttackList_[i])
                    )
                );
            }
        }
    }
    else
    {
        forAll(angleOfAttackList_, i)
        {
            if 
            (
                angleOfAttackList_[i] >= alphaDegStart
                and
                angleOfAttackList_[i] <= alphaDegStop
            )
            {
                newList.append
                (
                    liftCoefficientList_[i]*sin(degToRad(angleOfAttackList_[i]))
                  - dragCoefficientList_[i]*cos(degToRad(angleOfAttackList_[i]))
                );
            }
        }
    }
    return newList;
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
