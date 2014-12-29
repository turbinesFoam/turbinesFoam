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
    dict_.lookup("spanLength") >> spanLength_;
    dict_.lookup("spanDirection") >> spanDirection_;
    dict_.lookup("coefficientData") >> coefficientData_;
    dict_.lookup("freeStreamDirection") >> freeStreamDirection_;
    
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
        return 0.0;
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
    mesh_(mesh)
{
    velocity_ = vector::zero;
    forceVector_ = vector::zero;
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


void Foam::fv::actuatorLineElement::calculate
(
    vectorField Uin,
    volVectorField& forceField
)
{
    // Calculate local wind velocity
    scalar upstreamDistance = chordLength_;
    vector upstreamPoint = position_ - upstreamDistance*freeStreamDirection_;
    // Get velocity vector at upstream point
    vector inflowVelocity = Uin[mesh_.findCell(upstreamPoint)];
    
    // Calculate relative velocity (note these are not projected onto a
    // plane perpendicular to the chord and span direction)
    vector relativeVelocity = inflowVelocity - velocity_;
    
    // Calculate angle of attack (radians)
    scalar angleOfAttackRad = acos((-chordDirection_ & relativeVelocity)
                            / (mag(chordDirection_)*mag(relativeVelocity)));
    angleOfAttack_ = angleOfAttackRad/Foam::constant::mathematical::pi*180.0;
    
    // Lookup lift and drag coefficients
    liftCoefficient_ = 0.5;
    dragCoefficient_ = 0.05;
    
    // Calculate force per unit density
    scalar area = chordLength_*spanLength_;
    scalar magSqrU = magSqr(relativeVelocity);
    scalar lift = 0.5*area*liftCoefficient_*magSqrU;
    scalar drag = 0.5*area*dragCoefficient_*magSqrU;
    vector liftDirection = spanDirection_ ^ relativeVelocity;
    liftDirection /= mag(liftDirection);
    vector dragDirection = relativeVelocity/mag(relativeVelocity);
    forceVector_ = lift*liftDirection + drag*dragDirection;
    
    // Calculate force field 
    scalar epsilon = 3.5; // An estimate from SOWFA tutorial
    scalar projectionRadius = (epsilon*Foam::sqrt(Foam::log(1.0/0.001)));
    
    // Find the cells within the element's sphere of influence
    DynamicList<label> sphereCells;
    scalar sphereRadius = chordLength_ + projectionRadius;
    forAll(mesh_.cells(),cellI)
    {
        if (mag(mesh_.C()[cellI] - position_) <= sphereRadius)
        {
            sphereCells.append(cellI);
        }
    }
    
    if (debug)
    {
        Info<< "Calculating force contribution from actuatorLineElement " 
            << name_ << endl;
        Info<< "    inflowVelocity: " << inflowVelocity << endl;
        Info<< "    relativeVelocity: " << relativeVelocity << endl;
        Info<< "    angleOfAttack (degrees): " << angleOfAttack_ << endl;
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
    RM.xx() = Foam::sqr(axis.x()) + (1.0 - Foam::sqr(axis.x())) * Foam::cos(angle); 
    RM.xy() = axis.x() * axis.y() * (1.0 - Foam::cos(angle)) - axis.z() * Foam::sin(angle); 
    RM.xz() = axis.x() * axis.z() * (1.0 - Foam::cos(angle)) + axis.y() * Foam::sin(angle);
    RM.yx() = axis.x() * axis.y() * (1.0 - Foam::cos(angle)) + axis.z() * Foam::sin(angle); 
    RM.yy() = Foam::sqr(axis.y()) + (1.0 - Foam::sqr(axis.y())) * Foam::cos(angle);
    RM.yz() = axis.y() * axis.z() * (1.0 - Foam::cos(angle)) - axis.x() * Foam::sin(angle);
    RM.zx() = axis.x() * axis.z() * (1.0 - Foam::cos(angle)) - axis.y() * Foam::sin(angle);
    RM.zy() = axis.y() * axis.z() * (1.0 - Foam::cos(angle)) + axis.x() * Foam::sin(angle);
    RM.zz() = Foam::sqr(axis.z()) + (1.0 - Foam::sqr(axis.z())) * Foam::cos(angle);
    
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

    // Return the rotated point to its new location relative to the rotation point.
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


void Foam::fv::actuatorLineElement::pitch(scalar radians)
{
    rotate(position_, spanDirection_, radians, false);
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


Foam::vector& Foam::fv::actuatorLineElement::force()
{
    return forceVector_;
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
            name_ + ":force",
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

    const vectorField Uin(eqn.psi());
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

    const vectorField Uin(eqn.psi());
    calculate(Uin, forceI);

    // Add force to total actuator line force
    force += forceI;
}


// ************************************************************************* //
