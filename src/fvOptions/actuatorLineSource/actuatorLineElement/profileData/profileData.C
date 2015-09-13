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
                
        if (coefficientData[i].size() == 4)
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
    Re_(1),
    refRe_(1)
{
    read();
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


// ************************************************************************* //
