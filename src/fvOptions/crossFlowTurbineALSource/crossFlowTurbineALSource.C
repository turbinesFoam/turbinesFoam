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

#include "crossFlowTurbineALSource.H"
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
    defineTypeNameAndDebug(crossFlowTurbineALSource, 0);
    addToRunTimeSelectionTable
    (
        option,
        crossFlowTurbineALSource,
        dictionary
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fv::crossFlowTurbineALSource::rotateVector
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
    RM.xx() = Foam::sqr(axis.x()) + (1.0 - Foam::sqr(axis.x())) * Foam::cos(angle); 
    RM.xy() = axis.x() * axis.y() * (1.0 - Foam::cos(angle)) - axis.z() * Foam::sin(angle); 
    RM.xz() = axis.x() * axis.z() * (1.0 - Foam::cos(angle)) + axis.y() * Foam::sin(angle);
    RM.yx() = axis.x() * axis.y() * (1.0 - Foam::cos(angle)) + axis.z() * Foam::sin(angle); 
    RM.yy() = Foam::sqr(axis.y()) + (1.0 - Foam::sqr(axis.y())) * Foam::cos(angle);
    RM.yz() = axis.y() * axis.z() * (1.0 - Foam::cos(angle)) - axis.x() * Foam::sin(angle);
    RM.zx() = axis.x() * axis.z() * (1.0 - Foam::cos(angle)) - axis.y() * Foam::sin(angle);
    RM.zy() = axis.y() * axis.z() * (1.0 - Foam::cos(angle)) + axis.x() * Foam::sin(angle);
    RM.zz() = Foam::sqr(axis.z()) + (1.0 - Foam::sqr(axis.z())) * Foam::cos(angle);
    
    // Rotation matrices make a rotation about the origin, so need to subtract 
    // rotation point off the point to be rotated.
    vectorToRotate -= rotationPoint;

    // Perform the rotation.
    vectorToRotate = RM & vectorToRotate;

    // Return the rotated point to its new location relative to the rotation point.
    vectorToRotate += rotationPoint;
}


void Foam::fv::crossFlowTurbineALSource::createCoordinateSystem()
{
    // Construct the local rotor coordinate system
    freeStreamDirection_ = freeStreamVelocity_/mag(freeStreamVelocity_);
    radialDirection_ = axis_^freeStreamDirection_;
    radialDirection_ = radialDirection_/mag(radialDirection_);
}


Foam::tmp<Foam::vectorField> Foam::fv::crossFlowTurbineALSource::inflowVelocity
(
    const volVectorField& U
) const
{
    return U.internalField();
}


void Foam::fv::crossFlowTurbineALSource::createBlades()
{
    int nBlades = nBlades_;
    blades_.setSize(nBlades);
    int nElements;
    word profileName;
    List<List<scalar> > elementData;
    List<List<scalar> > profileData;
    word modelType = "actuatorLineSource";
    
    const dictionary& profilesSubDict(coeffs_.subDict("profiles"));
    
    for (int i = 0; i < nBlades_; i++)
    {
        word& bladeName = bladeNames_[i];
        // Create dictionary items for this blade
        dictionary bladeSubDict;
        bladeSubDict = bladesDict_.subDict(bladeName);
        bladeSubDict.lookup("nElements") >> nElements;
        bladeSubDict.lookup("profile") >> profileName;
        bladeSubDict.lookup("elementData") >> elementData;
        profilesSubDict.subDict(profileName).lookup("data") >> profileData;
        
        bladeSubDict.add("fieldNames", coeffs_.lookup("fieldNames"));
        bladeSubDict.add("coefficientData", profileData);
        bladeSubDict.add("tipEffect", tipEffect_);
        bladeSubDict.add("freeStreamVelocity", freeStreamVelocity_);
        
        if (debug)
        {
            Info<< "Creating actuator line blade " << bladeName << endl;
            Info<< "Blade has " << nElements << " elements" << endl;
            Info<< "Blade profile: " << profileName << endl;
            Info<< "Element data:" << endl;
            Info<< elementData << endl << endl;
            Info<< "Profile sectional coefficient data:" << endl;
            Info<< profileData << endl << endl;
        }
        
        // Convert element data into actuator line element geometry
        label nGeomPoints = elementData.size();
        List<List<List<scalar> > > elementGeometry(nGeomPoints);
        List<vector> initialVelocities(nGeomPoints, vector::one);
        for (int j = 0; j < nGeomPoints; j++)
        {
            // Read CFTAL dict data
            scalar axialDistance = elementData[j][0];
            scalar radius = elementData[j][1];
            scalar azimuthDegrees = elementData[j][2];
            scalar azimuthRadians = azimuthDegrees/180.0*mathematical::pi;
            scalar chordLength = elementData[j][3];
            scalar chordMount = elementData[j][4];
            
            // Set sizes for actuatorLineSource elementGeometry lists
            elementGeometry[j].setSize(4);
            elementGeometry[j][0].setSize(3);
            elementGeometry[j][1].setSize(3);
            elementGeometry[j][2].setSize(1);
            elementGeometry[j][3].setSize(1);
            
            // Create geometry point for AL source at origin
            vector point = origin_;
            // Move along axis
            point += axialDistance*axis_;
            scalar chordDisplacement = (0.5 - chordMount)*chordLength;
            point += chordDisplacement*freeStreamDirection_;
            point += radius*radialDirection_;
            initialVelocities[j] = -freeStreamDirection_*omega_*radius;
            // Rotate according to azimuth value
            rotateVector(point, origin_, axis_, azimuthRadians);
            rotateVector(initialVelocities[j], origin_, axis_, azimuthRadians);
            elementGeometry[j][0][0] = point.x(); // x location of geom point
            elementGeometry[j][0][1] = point.y(); // y location of geom point
            elementGeometry[j][0][2] = point.z(); // z location of geom point
            
            // Set span directions for AL source
            elementGeometry[j][1][0] = axis_.x(); // x component of span direction
            elementGeometry[j][1][1] = axis_.y(); // y component of span direction
            elementGeometry[j][1][2] = axis_.z(); // z component of span direction
            
            // Set chord length
            elementGeometry[j][2][0] = chordLength;
            
            // Set pitch
            scalar pitch = elementData[j][5];
            pitch += azimuthDegrees;
            elementGeometry[j][3][0] = pitch;
        }
        
        if (debug)
        {
            Info<< "Converted element geometry:" << endl << elementGeometry 
                << endl;
        }
        
        bladeSubDict.add("elementGeometry", elementGeometry);
        bladeSubDict.add("initialVelocities", initialVelocities);
        
        dictionary dict;
        dict.add("actuatorLineSourceCoeffs", bladeSubDict);
        dict.add("type", "actuatorLineSource");
        dict.add("active", dict_.lookup("active"));
        dict.add("selectionMode", dict_.lookup("selectionMode"));
        dict.add("cellSet", dict_.lookup("cellSet"));
        
        actuatorLineSource* blade = new actuatorLineSource
        (
            bladeName, 
            modelType, 
            dict, 
            mesh_
        );
        
        blades_.set(i, blade);
    }
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::fv::crossFlowTurbineALSource::crossFlowTurbineALSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    option(name, modelType, dict, mesh),
    time_(mesh.time()),
    rhoRef_(1.0),
    omega_(0.0),
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
    )
{
    read(dict);
    createCoordinateSystem();
    createBlades();
    lastRotationTime_ = time_.value();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::crossFlowTurbineALSource::~crossFlowTurbineALSource()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::vector& Foam::fv::crossFlowTurbineALSource::force()
{
    return force_;
}


Foam::volVectorField& Foam::fv::crossFlowTurbineALSource::forceField()
{
    return forceField_;
}

void Foam::fv::crossFlowTurbineALSource::rotate()
{
    scalar deltaT = time_.deltaT().value();
    scalar radians = omega_*deltaT;
    forAll(blades_, i)
    {
        blades_[i].rotate(origin_, axis_, radians);
    }
    
    if (debug)
    {
        Info<< "Rotating " << name_ << " " << radians << " radians" 
            << endl << endl;
    }
    
    lastRotationTime_ = time_.value();
}


void Foam::fv::crossFlowTurbineALSource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    // Rotate the turbine if time value has changed
    if (time_.value() != lastRotationTime_)
    {
        rotate();
    }

    // Zero out force vector and field
    forceField_ *= 0;
    force_ *= 0;
    
    // Read the reference density for incompressible flow
    //coeffs_.lookup("rhoRef") >> rhoRef_;

    // Add source for all actuator lines
    forAll(blades_, i)
    {
        blades_[i].addSup(eqn, fieldI);
        forceField_ += blades_[i].forceField();
        force_ += blades_[i].force();
    }
}


void Foam::fv::crossFlowTurbineALSource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    // Rotate the turbine if time value has changed
    if (time_.value() != lastRotationTime_)
    {
        rotate();
    }

    // Zero out force vector and field
    forceField_ *= 0;
    force_ *= 0;

    // Add source for all actuator lines
    forAll(blades_, i)
    {
        blades_[i].addSup(rho, eqn, fieldI);
        forceField_ += blades_[i].forceField();
        force_ += blades_[i].force();
    }
}


void Foam::fv::crossFlowTurbineALSource::writeData(Ostream& os) const
{
    os << indent << name_ << endl;
    dict_.write(os);
}


void Foam::fv::crossFlowTurbineALSource::printCoeffs() const
{
    Info<< "Number of blades: " << nBlades_ << endl;
}


bool Foam::fv::crossFlowTurbineALSource::read(const dictionary& dict)
{
    if (option::read(dict))
    {
        
        coeffs_.lookup("fieldNames") >> fieldNames_;
        applied_.setSize(fieldNames_.size(), false);

        // Read coordinate system/geometry invariant properties
        coeffs_.lookup("origin") >> origin_;
        coeffs_.lookup("axis") >> axis_;
        coeffs_.lookup("freeStreamVelocity") >> freeStreamVelocity_;
        coeffs_.lookup("tipSpeedRatio") >> tipSpeedRatio_;
        coeffs_.lookup("rotorRadius") >> rotorRadius_;
        omega_ = tipSpeedRatio_*mag(freeStreamVelocity_)/rotorRadius_;

        // Get blade information
        bladesDict_ = coeffs_.subDict("blades");
        nBlades_ = bladesDict_.keys().size();
        bladeNames_ = bladesDict_.toc();

        coeffs_.lookup("tipEffect") >> tipEffect_;
        
        if (debug)
        {
            Info<< "Debugging on" << endl;
            Info<< "Cross-flow turbine properties:" << endl;
            printCoeffs();
        }

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
