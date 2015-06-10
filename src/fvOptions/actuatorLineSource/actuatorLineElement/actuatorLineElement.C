/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "actuatorLineElement.H"
#include "addToRunTimeSelectionTable.H"
#include "geometricOneField.H"
#include "fvMatrices.H"
#include "syncTools.H"
#include "interpolationCellPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(actuatorLineElement, 0);
    defineRunTimeSelectionTable(actuatorLineElement, dictionary);
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::fv::actuatorLineElement::read()
{
    // Parse dictionary
    dict_.lookup("position") >> position_;
    dict_.lookup("chordLength") >> chordLength_;
    dict_.lookup("chordDirection") >> chordDirection_;
    dict_.lookup("chordMount") >> chordMount_;
    dict_.lookup("spanLength") >> spanLength_;
    dict_.lookup("spanDirection") >> spanDirection_;
    dict_.lookup("coefficientData") >> coefficientData_;
    dict_.lookup("freeStreamVelocity") >> freeStreamVelocity_;
    freeStreamDirection_ = freeStreamVelocity_/mag(freeStreamVelocity_);
    
    // Create lists from coefficient data
    angleOfAttackList_.setSize(coefficientData_.size());
    liftCoefficientList_.setSize(coefficientData_.size());
    dragCoefficientList_.setSize(coefficientData_.size());
    forAll(coefficientData_, i)
    {
        angleOfAttackList_[i] = coefficientData_[i][0];
        liftCoefficientList_[i] = coefficientData_[i][1];
        dragCoefficientList_[i] = coefficientData_[i][2];
    }
    
    // Create dynamic stall model if found
    if (dict_.found("dynamicStall"))
    {
        dictionary dsDict = dict_.subDict("dynamicStall");
        word dsName = dsDict.lookup("dynamicStallModel");
        dynamicStall_ = dynamicStallModel::New
        (
            dsDict, 
            dsName, 
            mesh_.time()
        );
        dsDict.lookup("active") >> dynamicStallActive_;
    }
    
    // Read nu from object registry
    const dictionary& transportProperties = mesh_.lookupObject<IOdictionary>
    (
        "transportProperties"
    );
    dimensionedScalar nu;
    transportProperties.lookup("nu") >> nu;
    nu_ = nu.value();
    
    if (debug)
    {
       Info<< "actuatorLineElement properties:" << endl;
       Info<< "Position: " << position_ << endl;
       Info<< "chordLength: " << chordLength_ << endl;
       Info<< "chordDirection: " << chordDirection_ << endl;
       Info<< "spanLength: " << spanLength_ << endl;
       Info<< "spanDirection: " << spanDirection_ << endl;
       Info<< "coefficientData: " << coefficientData_ << endl;
    }
}


Foam::scalar Foam::fv::actuatorLineElement::interpolate
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
        if (debug)
        {
            Info<< "    Interpolation failed" << endl;
        }
        return 0.0;
    }
}


void Foam::fv::actuatorLineElement::lookupCoefficients()
{
    liftCoefficient_ = interpolate
    (
        angleOfAttack_, 
        angleOfAttackList_,
        liftCoefficientList_
    );
    dragCoefficient_ = interpolate
    (
        angleOfAttack_, 
        angleOfAttackList_,
        dragCoefficientList_
    );
}


Foam::scalar Foam::fv::actuatorLineElement::calcProjectionEpsilon()
{
    const scalarField& V = mesh_.V();
    label posCellI = mesh_.findCell(position_);
    // Projection width based on local cell size (from Troldborg (2008))
    scalar epsilon = 2*Foam::cbrt(V[posCellI]);
    
    if (epsilon > (chordLength_/2.0))
    {
        return epsilon;
    }
    else
    {
        return chordLength_/2.0;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::actuatorLineElement::actuatorLineElement
(
    const word& name,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    dict_(dict),
    name_(name),
    mesh_(mesh),
    velocity_(vector::zero),
    forceVector_(vector::zero),
    relativeVelocity_(vector::zero),
    angleOfAttack_(0.0),
    dynamicStallActive_(false),
    omega_(0.0),
    chordMount_(0.25)
{
    read();
}

// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

Foam::fv::actuatorLineElement::~actuatorLineElement()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::word& Foam::fv::actuatorLineElement::name() const
{
    return name_;
}


Foam::vector& Foam::fv::actuatorLineElement::position()
{
    return position_;
}


Foam::vector& Foam::fv::actuatorLineElement::relativeVelocity()
{
    return relativeVelocity_;
}


Foam::scalar& Foam::fv::actuatorLineElement::angleOfAttack()
{
    return angleOfAttack_;
}


Foam::scalar& Foam::fv::actuatorLineElement::angleOfAttackGeom()
{
    return angleOfAttackGeom_;
}


Foam::scalar& Foam::fv::actuatorLineElement::liftCoefficient()
{
    return liftCoefficient_;
}


Foam::scalar& Foam::fv::actuatorLineElement::dragCoefficient()
{
    return dragCoefficient_;
}


void Foam::fv::actuatorLineElement::calculate
(
    const volVectorField& Uin,
    volVectorField& forceField
)
{
    scalar pi = Foam::constant::mathematical::pi;
    if (debug)
    {
        Info<< "Calculating force contribution from actuatorLineElement " 
            << name_ << endl;
        Info<< "    position: " << position_ << endl;
        Info<< "    chordDirection: " << chordDirection_ << endl;
        Info<< "    spanDirection: " << spanDirection_ << endl;
        Info<< "    elementVelocity: " << velocity_ << endl;
    }
    
    // Find local flow velocity upstream
    scalar upstreamDistance = chordLength_*0.0;
    vector upstreamPoint = position_ - upstreamDistance*freeStreamDirection_;
    interpolationCellPoint<vector> UInterp(Uin);
    label upstreamCellI = mesh_.findCell(upstreamPoint);
    vector inflowVelocity = UInterp.interpolate
    (
        upstreamPoint, 
        upstreamCellI
    );
    
    if (debug)
    {
        Info<< "    inflowVelocity: " << inflowVelocity << endl;
    }
    
    // Calculate relative velocity (note these are not projected onto a
    // plane perpendicular to the chord and span direction)
    relativeVelocity_ = inflowVelocity - velocity_;
    Re_ = mag(relativeVelocity_)*chordLength_/nu_;
    
    if (debug)
    {
        Info<< "    Relative upstream point: " << (upstreamPoint - position_)
            << endl;
        Info<< "    relativeVelocity: " << relativeVelocity_ << endl;
        Info<< "    Reynolds number: " << Re_ << endl;
    }
    
    // Calculate vector normal to chord--span plane
    vector planformNormal = -chordDirection_ ^ spanDirection_;
    planformNormal /= mag(planformNormal);
    
    // Calculate angle of attack (radians)
    scalar angleOfAttackRad = asin((planformNormal & relativeVelocity_)
                            / (mag(planformNormal)*mag(relativeVelocity_)));
    vector relVelGeom = freeStreamVelocity_ - velocity_;
    angleOfAttackGeom_ = asin((planformNormal & relVelGeom)
                       / (mag(planformNormal)*mag(relVelGeom)));
    angleOfAttackGeom_ *= 180.0/pi;
                            
    if (debug)
    {
        Info<< "    Angle of attack (uncorrected, degrees): " 
            << angleOfAttackRad/pi*180.0 << endl;
        Info<< "    Geometric angle of attack (degrees): "
            << angleOfAttackGeom_ << endl;
    }
    
    // Apply flow curvature correction to angle of attack
    angleOfAttackRad += omega_*(chordMount_ - 0.25)
                      * chordLength_/mag(relativeVelocity_);
    angleOfAttackRad += omega_*chordLength_/(4*mag(relativeVelocity_));
    
    // Calculate angle of attack in degrees
    angleOfAttack_ = angleOfAttackRad/pi*180.0;
    
    if (debug)
    {
        Info<< "    Angle of attack (corrected, degrees): " 
            << angleOfAttack_ << endl;
    }
    
    // Lookup lift and drag coefficients
    lookupCoefficients();
    
    // Correct coefficients with dynamic stall model
    if (dynamicStallActive_)
    {
        dynamicStall_->correct
        (
            mag(relativeVelocity_),
            angleOfAttack_,
            liftCoefficient_,
            dragCoefficient_,
            angleOfAttackList_,
            liftCoefficientList_,
            dragCoefficientList_
        );
    }
    
    // Calculate force per unit density
    scalar area = chordLength_*spanLength_;
    scalar magSqrU = magSqr(relativeVelocity_);
    scalar lift = 0.5*area*liftCoefficient_*magSqrU;
    scalar drag = 0.5*area*dragCoefficient_*magSqrU;
    vector liftDirection = relativeVelocity_ ^ spanDirection_;
    liftDirection /= mag(liftDirection);
    vector dragDirection = relativeVelocity_/mag(relativeVelocity_);
    forceVector_ = lift*liftDirection + drag*dragDirection;
    
    // Calculate projection width
    scalar epsilon = calcProjectionEpsilon();
    scalar projectionRadius = (epsilon*Foam::sqrt(Foam::log(1.0/0.001)));
    
    // Apply force to the cells within the element's sphere of influence
    scalar sphereRadius = chordLength_ + projectionRadius;
    forAll(mesh_.cells(), cellI)
    {
        scalar dis = mag(mesh_.C()[cellI] - position_);
        if (dis <= sphereRadius)
        {
            scalar factor = Foam::exp(-Foam::sqr(dis/epsilon))
                          / (Foam::pow(epsilon, 3)
                          * Foam::pow(Foam::constant::mathematical::pi, 1.5));
            // forceField is opposite forceVector
            forceField[cellI] += -forceVector_*factor;
        }
    }
    
    
    if (debug)
    {
        Info<< "    epsilon: " << epsilon << endl;
        Info<< "    sphereRadius: " << sphereRadius << endl;
        Info<< "    force (per unit density): " << forceVector_ << endl 
            << endl;
    }
}


void Foam::fv::actuatorLineElement::rotate
(
    vector rotationPoint, 
    vector axis, 
    scalar radians,
    bool rotateVelocity=true
)
{
    // Declare and define the rotation matrix (from SOWFA)
    tensor RM;
    scalar angle = radians;
    RM.xx() = Foam::sqr(axis.x()) 
            + (1.0 - Foam::sqr(axis.x())) * Foam::cos(angle); 
    RM.xy() = axis.x() * axis.y() 
            * (1.0 - Foam::cos(angle)) - axis.z() * Foam::sin(angle); 
    RM.xz() = axis.x() * axis.z() 
            * (1.0 - Foam::cos(angle)) + axis.y() * Foam::sin(angle);
    RM.yx() = axis.x() * axis.y() 
            * (1.0 - Foam::cos(angle)) + axis.z() * Foam::sin(angle); 
    RM.yy() = Foam::sqr(axis.y()) 
            + (1.0 - Foam::sqr(axis.y())) * Foam::cos(angle);
    RM.yz() = axis.y() * axis.z() 
            * (1.0 - Foam::cos(angle)) - axis.x() * Foam::sin(angle);
    RM.zx() = axis.x() * axis.z() 
            * (1.0 - Foam::cos(angle)) - axis.y() * Foam::sin(angle);
    RM.zy() = axis.y() * axis.z() 
            * (1.0 - Foam::cos(angle)) + axis.x() * Foam::sin(angle);
    RM.zz() = Foam::sqr(axis.z()) 
            + (1.0 - Foam::sqr(axis.z())) * Foam::cos(angle);
    
    if (debug)
    {
        Info<< "Rotating actuatorLineElement: " << name_ << endl;
        Info<< "Rotation point: " << rotationPoint << endl;
        Info<< "Rotation axis: " << axis << endl;
        Info<< "Rotation angle (radians): " << radians << endl;
        Info<< "Rotation matrix:" << endl << RM << endl;
        Info<< "Initial position: " << position_ << endl;
        Info<< "Initial chordDirection: " << chordDirection_ << endl;
        Info<< "Initial spanDirection: " << spanDirection_ << endl;
        Info<< "Initial velocity: " << velocity_ << endl;
    }

    // Rotation matrices make a rotation about the origin, so need to subtract 
    // rotation point off the point to be rotated.
    vector point = position_;
    point -= rotationPoint;

    // Perform the rotation.
    point = RM & point;

    // Return the rotated point to its new location relative to the rotation 
    // point
    point += rotationPoint;
    
    // Set the position of the element
    position_ = point;
    
    // Rotate the span and chord vectors of the element
    chordDirection_ = RM & chordDirection_;
    spanDirection_ = RM & spanDirection_;
    
    // Rotate the element's velocity vector if specified
    if (rotateVelocity) velocity_ = RM & velocity_;
    
    if (debug)
    {
        Info<< "Final position: " << position_ << endl;
        Info<< "Final chordDirection: " << chordDirection_ << endl;
        Info<< "Final spanDirection: " << spanDirection_ << endl;
        Info<< "Final velocity: " << velocity_ << endl << endl;
    }
}


void Foam::fv::actuatorLineElement::pitch
(
    scalar radians, 
    scalar chordFraction
)
{
    vector rotationPoint = position_;
    rotationPoint += chordDirection_*(chordMount_ - chordFraction);
    rotate(rotationPoint, spanDirection_, radians, false);
}


void Foam::fv::actuatorLineElement::translate(vector translationVector)
{
    position_ += translationVector;
}


void Foam::fv::actuatorLineElement::setVelocity(vector velocity)
{
    if (debug) 
    {
        Info<< "Changing velocity of " << name_ << " from " 
            << velocity_ << " to " << velocity << endl << endl;
    }
    velocity_ = velocity;
}


void Foam::fv::actuatorLineElement::setSpeed(scalar speed)
{
    if (mag(velocity_) > 0)
    {
        velocity_ /= mag(velocity_);
        velocity_ *= speed;
    }
}


void Foam::fv::actuatorLineElement::setSpeed
(
    vector point,
    vector axis,
    scalar omega
)
{
    if (debug)
    {
        Info<< "Setting speed of " << name_ << " from rotation" << endl;
        Info<< "    Initial velocity: " << velocity_ << endl;
    }
    
    // First find radius from axis to element position -- formula from 
    // http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
    vector point2 = point + axis;
    scalar radius = mag((position_ - point) ^ (position_ - point2))
                  / mag(point2 - point);
    scalar speed = omega*radius;
    setSpeed(speed);
    // Also set omega for flow curvature correction
    setOmega(omega);
    
    if (debug)
    {
        Info<< "    Radius: " << radius << endl;
        Info<< "    Final velocity: " << velocity_ << endl;
    }
}


void Foam::fv::actuatorLineElement::scaleVelocity(scalar scale)
{
    velocity_ *= scale;
}


Foam::vector& Foam::fv::actuatorLineElement::force()
{
    return forceVector_;
}


Foam::vector Foam::fv::actuatorLineElement::moment(vector point)
{
    // Calculate radius vector
    vector radius = position_ - point;
    vector moment = radius ^ forceVector_;
    return moment;
}


void Foam::fv::actuatorLineElement::addSup
(
    fvMatrix<vector>& eqn,
    volVectorField& force
)
{
    volVectorField forceI
    (
        IOobject
        (
            "force." + name_,
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedVector
        (
            "zero",
            eqn.dimensions()/dimVolume,
            vector::zero
        )
    );

    // Read the reference density for incompressible flow
    //coeffs_.lookup("rhoRef") >> rhoRef_;

    const volVectorField& Uin(eqn.psi());
    calculate(Uin, forceI);

    // Add source to rhs of eqn
    force += forceI;
}


void Foam::fv::actuatorLineElement::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    volVectorField& force
)
{
    volVectorField forceI
    (
        IOobject
        (
            name_ + "Force",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedVector
        (
            "zero",
            eqn.dimensions()/dimVolume,
            vector::zero
        )
    );

    const volVectorField& Uin(eqn.psi());
    calculate(Uin, forceI);

    // Add force to total actuator line force
    force += forceI;
}


void Foam::fv::actuatorLineElement::addTurbulence
(   
    fvMatrix<scalar>& eqn,
    word fieldName
)
{
    volScalarField turbulence
    (
        IOobject
        (
            "turbulence." + name_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar
        (
            "zero", 
            eqn.dimensions()/dimVolume, 
            0.0
        )
    );
    
    // Calculate projection radius
    scalar epsilon = calcProjectionEpsilon();
    scalar projectionRadius = (epsilon*Foam::sqrt(Foam::log(1.0/0.001)));
    
    // Calculate TKE injection rate
    scalar k = 0.1*mag(dragCoefficient_);
    
    // Add turbulence to the cells within the element's sphere of influence
    scalar sphereRadius = chordLength_ + projectionRadius;
    forAll(mesh_.cells(), cellI)
    {
        scalar dis = mag(mesh_.C()[cellI] - position_);
        if (dis <= sphereRadius)
        {
            scalar factor = Foam::exp(-Foam::sqr(dis/epsilon))
                          / (Foam::pow(epsilon, 3)
                          * Foam::pow(Foam::constant::mathematical::pi, 1.5));
            if (fieldName == "k")
            {
                turbulence[cellI] = factor*k;
            }
            else if (fieldName == "epsilon")
            {
                turbulence[cellI] = factor*Foam::pow(k, 1.5)
                                  * 0.09/(chordLength_/10.0);
            }
        }
    }
    
    eqn += turbulence;
}


void Foam::fv::actuatorLineElement::setDynamicStallActive(bool active)
{
    dynamicStallActive_ = active;
}


void Foam::fv::actuatorLineElement::setOmega(scalar omega)
{
    omega_ = omega;
}

// ************************************************************************* //
