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

#include "profileData.H"
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


void Foam::profileData::read()
{
    // Read reference Reynolds number, and if present turn on Reynolds number
    // corrections
    ReRef_ = dict_.lookupOrDefault("Re", VSMALL);
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


void Foam::profileData::analyze()
{
    calcStaticStallAngle();
    calcZeroLiftDragCoeff();
    calcZeroLiftAngleOfAttack();
    calcNormalCoeffSlope();
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::profileData::profileData
(
    const word& name,
    const dictionary& dict
)
:
    name_(name),
    dict_(dict),
    Re_(VSMALL),
    ReRef_(VSMALL),
    correctRe_(false),
    staticStallAngle_(0.0),
    zeroLiftDragCoeff_(0.0),
    zeroLiftAngleOfAttack_(0.0),
    normalCoeffSlope_(0.0)
{
    read();
    // Analyze if the data is not from a bluff body
    if (Foam::sum(Foam::mag(liftCoefficientList_)) > VSMALL)
    {
        analyze();
    }
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::profileData>
Foam::profileData::New()
{
    return autoPtr<profileData>(new profileData);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::profileData::~profileData()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

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


void Foam::profileData::updateRe(scalar Re)
{
    if (correctRe_ and Re != Re_)
    {
        // Correct drag coefficients
        scalar K = pow((Re/ReRef_), -0.5);
        dragCoefficientList_ = K*dragCoefficientListOrg_;
        
        // Correct lift coefficients
        scalar n = 0.125;
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
        
        // Recalculate static stall angle, etc.
        analyze();
    }
    
    Re_ = Re;
}


Foam::List<scalar>& Foam::profileData::angleOfAttackList()
{
    return angleOfAttackList_;
}


Foam::List<scalar>& Foam::profileData::liftCoefficientList()
{
    return liftCoefficientList_;
}


Foam::List<scalar>& Foam::profileData::dragCoefficientList()
{
    return dragCoefficientList_;
}


Foam::List<scalar>& Foam::profileData::momentCoefficientList()
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
    return newList;
}


Foam::List<scalar> Foam::profileData::liftCoefficientList
(
    scalar alphaDegStart,
    scalar alphaDegStop
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
            newList.append(liftCoefficientList_[i]);
        }
    }
    return newList;
}


Foam::List<scalar> Foam::profileData::dragCoefficientList
(
    scalar alphaDegStart,
    scalar alphaDegStop
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
            newList.append(dragCoefficientList_[i]);
        }
    }
    return newList;
}


Foam::List<scalar> Foam::profileData::momentCoefficientList
(
    scalar alphaDegStart,
    scalar alphaDegStop
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
            newList.append(momentCoefficientList_[i]);
        }
    }
    return newList;
}


Foam::List<scalar> Foam::profileData::normalCoefficientList
(
    scalar alphaDegStart,
    scalar alphaDegStop
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
            newList.append
            (
                liftCoefficientList_[i]*cos(degToRad(angleOfAttackList_[i]))
              - dragCoefficientList_[i]*sin(degToRad(angleOfAttackList_[i]))
            );
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
    return newList;
}


Foam::scalar Foam::profileData::staticStallAngleRad()
{
    return degToRad(staticStallAngle_);
}


Foam::scalar Foam::profileData::zeroLiftDragCoeff()
{
    return zeroLiftDragCoeff_;
}


Foam::scalar Foam::profileData::zeroLiftAngleOfAttack()
{
    return zeroLiftAngleOfAttack_;
}


Foam::scalar Foam::profileData::normalCoeffSlope()
{
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
