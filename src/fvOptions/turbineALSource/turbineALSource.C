/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

#include "turbineALSource.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMatrices.H"
#include "geometricOneField.H"
#include "syncTools.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(turbineALSource, 0);
    addToRunTimeSelectionTable
    (
        option,
        turbineALSource,
        dictionary
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fv::turbineALSource::rotateVector
(
    vector& vectorToRotate,
    vector rotationPoint, 
    vector axis,
    scalar radians
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
    
    // Rotation matrices make a rotation about the origin, so need to subtract 
    // rotation point off the point to be rotated.
    vectorToRotate -= rotationPoint;

    // Perform the rotation.
    vectorToRotate = RM & vectorToRotate;

    // Return the rotated point to its new location relative to the rotation 
    // point
    vectorToRotate += rotationPoint;
}


void Foam::fv::turbineALSource::createCoordinateSystem()
{
    // Should be unique for each type of turbine
}


void Foam::fv::turbineALSource::createBlades()
{
    // Should be unique for each type of turbine
}


void Foam::fv::turbineALSource::createOutputFile()
{
    fileName dir;
    
    if (Pstream::parRun())
    {
        dir = time_.path()/"../postProcessing/turbines"
            / time_.timeName();
    }
    else
    {
        dir = time_.path()/"postProcessing/turbines"
            / time_.timeName();
    }
    
    if (!isDir(dir))
    {
        mkDir(dir);
    }

    outputFile_ = new OFstream(dir/name_ + ".csv");
    
    *outputFile_<< "time,angle_deg,tsr,cp,cd,ct" << endl;
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::fv::turbineALSource::turbineALSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    option(name, modelType, dict, mesh),
    time_(mesh.time()),
    lastRotationTime_(time_.value()),
    rhoRef_(1.0),
    omega_(0.0),
    angleDeg_(0.0),
    nBlades_(0),
    freeStreamVelocity_(vector::zero),
    tipEffect_(1.0),
    forceField_
    (
        IOobject
        (
            "force." + name_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector
        (
            "force", 
            dimForce/dimVolume/dimDensity, 
            vector::zero
        )
    ),
    frontalArea_(0.0),
    powerCoefficient_(0.0),
    dragCoefficient_(0.0),
    torqueCoefficient_(0.0)
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::turbineALSource::~turbineALSource()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::vector& Foam::fv::turbineALSource::force()
{
    return force_;
}


Foam::volVectorField& Foam::fv::turbineALSource::forceField()
{
    return forceField_;
}


Foam::scalar& Foam::fv::turbineALSource::torque()
{
    return torque_;
}


void Foam::fv::turbineALSource::rotate()
{
}


void Foam::fv::turbineALSource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    // Should be unique for each turbine type
}


void Foam::fv::turbineALSource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    // Should be unique for each turbine type
}


void Foam::fv::turbineALSource::addSup
(
    fvMatrix<scalar>& eqn,
    const label fieldI
)
{
    // Should be unique for each turbine type
}


void Foam::fv::turbineALSource::printCoeffs() const
{
    Info<< "Number of blades: " << nBlades_ << endl;
}


void Foam::fv::turbineALSource::writePerf()
{
    *outputFile_<< time_.value() << "," << angleDeg_ << "," 
                << tipSpeedRatio_ << "," << powerCoefficient_ << "," 
                << dragCoefficient_ << "," << torqueCoefficient_ << endl;
}


void Foam::fv::turbineALSource::writeData(Ostream& os) const
{
    os  << indent << name_ << endl;
    dict_.write(os);
}


bool Foam::fv::turbineALSource::read(const dictionary& dict)
{
    if (option::read(dict))
    {
        return true;
    }
    else
    {
        return false;
    }
}

// ************************************************************************* //
