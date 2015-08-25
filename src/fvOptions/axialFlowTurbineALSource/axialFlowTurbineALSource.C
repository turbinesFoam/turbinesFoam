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

#include "axialFlowTurbineALSource.H"
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
    defineTypeNameAndDebug(axialFlowTurbineALSource, 0);
    addToRunTimeSelectionTable
    (
        option,
        axialFlowTurbineALSource,
        dictionary
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fv::axialFlowTurbineALSource::createCoordinateSystem()
{
    // Make sure axis is a unit vector
    axis_ /= mag(axis_);
    
    // Free stream direction is a unit vector
    freeStreamDirection_ = freeStreamVelocity_/mag(freeStreamVelocity_);
    
    // Radial direction is vertical direction
    verticalDirection_ /= mag(verticalDirection_);
    radialDirection_ = verticalDirection_;
    
    // Calculate initial azimuthal or tangential direction
    azimuthalDirection_ = axis_ ^ verticalDirection_;
    azimuthalDirection_ /= mag(azimuthalDirection_);
}


void Foam::fv::axialFlowTurbineALSource::createBlades()
{
    int nBlades = nBlades_;
    blades_.setSize(nBlades);
    int nElements;
    List<List<scalar> > elementData;
    word modelType = "actuatorLineSource";
    List<scalar> frontalAreas(nBlades); // frontal area from each blade
    
    for (int i = 0; i < nBlades_; i++)
    {
        word& bladeName = bladeNames_[i];
        // Create dictionary items for this blade
        dictionary bladeSubDict = bladesDict_.subDict(bladeName);
        bladeSubDict.lookup("nElements") >> nElements;
        bladeSubDict.lookup("elementData") >> elementData;
        scalar azimuthalOffset = bladeSubDict.lookupOrDefault
        (
            "azimuthalOffset",
            0.0
        );
        
        bladeSubDict.add("freeStreamVelocity", freeStreamVelocity_);
        bladeSubDict.add("fieldNames", coeffs_.lookup("fieldNames"));
        bladeSubDict.add("tipEffect", tipEffect_);
        bladeSubDict.add("profileData", profileData_);
        
        if (debug)
        {
            Info<< "Creating actuator line blade " << bladeName << endl;
            Info<< "Blade has " << nElements << " elements" << endl;
            Info<< "Element data:" << endl;
            Info<< elementData << endl << endl;
        }
        
        // Convert element data into actuator line element geometry
        label nGeomPoints = elementData.size();
        List<List<List<scalar> > > elementGeometry(nGeomPoints);
        List<vector> initialVelocities(nGeomPoints, vector::zero);
        // Frontal area for this blade
        scalar frontalArea = 0.0;
        scalar maxRadius = 0.0;
        forAll(elementData, j)
        {
            // Read AFTAL dict element data
            scalar axialDistance = elementData[j][0];
            scalar radius = elementData[j][1];
            scalar azimuthDegrees = elementData[j][2] + azimuthalOffset;
            scalar azimuthRadians = degToRad(azimuthDegrees);
            scalar chordLength = elementData[j][3];
            scalar chordMount = elementData[j][4];
            scalar pitch = elementData[j][5];
            
            // Find max radius for calculating frontal area
            if (radius > maxRadius)
            {
                maxRadius = radius;
            }
            
            // Set sizes for actuatorLineSource elementGeometry lists
            elementGeometry[j].setSize(6);
            elementGeometry[j][0].setSize(3);
            elementGeometry[j][1].setSize(3);
            elementGeometry[j][2].setSize(1);
            elementGeometry[j][3].setSize(3);
            elementGeometry[j][4].setSize(1);
            elementGeometry[j][5].setSize(1);
            
            // Create geometry point for AL source at origin
            vector point = origin_;
            // Move point along axial direction
            point += axialDistance*axis_;
            // Move along radial direction
            point += radius*radialDirection_;
            // Move along chord according to chordMount
            scalar chordDisplacement = (chordMount - 0.25)*chordLength;
            point -= chordDisplacement*azimuthalDirection_;
            // Set initial velocity of quarter chord
            scalar radiusCorr = sqrt(magSqr(chordMount - 0.25)*chordLength
                                     + magSqr(radius));
            vector initialVelocity = azimuthalDirection_*omega_*radiusCorr;
            scalar velAngle = atan2(((chordMount - 0.25)*chordLength), radius);
            rotateVector(initialVelocity, vector::zero, axis_, velAngle);
            initialVelocities[j] = initialVelocity;
            // Rotate point and initial velocity according to azimuth value
            rotateVector(point, origin_, axis_, azimuthRadians);
            rotateVector
            (
                initialVelocities[j],
                vector::zero,
                axis_,
                azimuthRadians
            );
            
            // Set point coordinates for AL source
            elementGeometry[j][0][0] = point.x(); // x location of geom point
            elementGeometry[j][0][1] = point.y(); // y location of geom point
            elementGeometry[j][0][2] = point.z(); // z location of geom point
            
            // Set span directions for AL source
            scalar spanSign = axis_ & freeStreamDirection_;
            vector spanDirection = spanSign*verticalDirection_;
            rotateVector(spanDirection, vector::zero, axis_, azimuthRadians);
            elementGeometry[j][1][0] = spanDirection.x();
            elementGeometry[j][1][1] = spanDirection.y();
            elementGeometry[j][1][2] = spanDirection.z();
            
            // Set chord length
            elementGeometry[j][2][0] = chordLength;
            
            // Set chord reference direction
            vector chordDirection = azimuthalDirection_;
            rotateVector(chordDirection, vector::zero, axis_, azimuthRadians);
            elementGeometry[j][3][0] = chordDirection.x(); 
            elementGeometry[j][3][1] = chordDirection.y(); 
            elementGeometry[j][3][2] = chordDirection.z();
            
            // Set chord mount
            elementGeometry[j][4][0] = chordMount;
            
            // Set pitch
            elementGeometry[j][5][0] = -pitch;
        }
        
        // Add frontal area to list
        frontalArea = mathematical::pi*magSqr(maxRadius);
        frontalAreas[i] = frontalArea;
        
        if (debug)
        {
            Info<< "Converted element geometry:" << endl << elementGeometry 
                << endl;
            Info<< "Frontal area from " << bladeName << ": " << frontalArea
                << endl;
        }
        
        bladeSubDict.add("elementGeometry", elementGeometry);
        bladeSubDict.add("initialVelocities", initialVelocities);
        bladeSubDict.add("dynamicStall", dynamicStallDict_);
        
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
    
    // Frontal area is calculated using defined rotorRadius rather than
    // detected from elementData
    frontalArea_ = mathematical::pi*magSqr(rotorRadius_);
    Info<< "Frontal area of " << name_ << ": " << frontalArea_ << endl;
}


void Foam::fv::axialFlowTurbineALSource::createHub()
{
    int nElements;
    List<List<scalar> > elementData;
    dictionary hubSubDict = hubDict_;
    
    hubDict_.lookup("nElements") >> nElements;
    hubDict_.lookup("elementData") >> elementData;
    
    // Convert element data into actuator line element geometry
    label nGeomPoints = elementData.size();
    List<List<List<scalar> > > elementGeometry(nGeomPoints);
    List<vector> initialVelocities(nGeomPoints, vector::zero);

    forAll(elementData, j)
    {
        // Read hub element data
        scalar axialDistance = elementData[j][0];
        scalar height = elementData[j][1];
        scalar diameter = elementData[j][2];
        
        // Set sizes for actuatorLineSource elementGeometry lists
        elementGeometry[j].setSize(6);
        elementGeometry[j][0].setSize(3);
        elementGeometry[j][1].setSize(3);
        elementGeometry[j][2].setSize(1);
        elementGeometry[j][3].setSize(3);
        elementGeometry[j][4].setSize(1);
        elementGeometry[j][5].setSize(1);
        
        // Create geometry point for AL source at origin
        vector point = origin_;
        // Move along axis
        point += axialDistance*axis_;
        // Move along vertical direction
        point += height*verticalDirection_;
        
        elementGeometry[j][0][0] = point.x(); // x location of geom point
        elementGeometry[j][0][1] = point.y(); // y location of geom point
        elementGeometry[j][0][2] = point.z(); // z location of geom point
        
        // Set span directions
        elementGeometry[j][1][0] = verticalDirection_.x();
        elementGeometry[j][1][1] = verticalDirection_.y();
        elementGeometry[j][1][2] = verticalDirection_.z();
        
        // Set chord length
        elementGeometry[j][2][0] = diameter;
        
        // Set chord reference direction
        elementGeometry[j][3][0] = freeStreamDirection_.x();
        elementGeometry[j][3][1] = freeStreamDirection_.y();
        elementGeometry[j][3][2] = freeStreamDirection_.z();
        
        // Set chord mount
        elementGeometry[j][4][0] = 0.25;
        
        // Set pitch
        elementGeometry[j][5][0] = 0.0;
    }
    
    hubSubDict.add("elementGeometry", elementGeometry);
    hubSubDict.add("initialVelocities", initialVelocities);
    hubSubDict.add("fieldNames", coeffs_.lookup("fieldNames"));
    hubSubDict.add("profileData", profileData_);
    hubSubDict.add("tipEffect", tipEffect_);
    hubSubDict.add("freeStreamVelocity", freeStreamVelocity_);
        
    dictionary dict;
    dict.add("actuatorLineSourceCoeffs", hubSubDict);
    dict.add("type", "actuatorLineSource");
    dict.add("active", dict_.lookup("active"));
    dict.add("selectionMode", dict_.lookup("selectionMode"));
    dict.add("cellSet", dict_.lookup("cellSet"));
    
    actuatorLineSource* hub = new actuatorLineSource
    (
        "hub", 
        "actuatorLineSource", 
        dict, 
        mesh_
    );
    
    hub_.set(hub);
}


void Foam::fv::axialFlowTurbineALSource::createTower()
{
    vector towerAxis = verticalDirection_;
    List<List<scalar> > elementData;
    dictionary towerSubDict = towerDict_;
    
    towerDict_.lookup("elementData") >> elementData;
    
    // Convert element data into actuator line element geometry
    label nGeomPoints = elementData.size();
    List<List<List<scalar> > > elementGeometry(nGeomPoints);
    List<vector> initialVelocities(nGeomPoints, vector::zero);

    forAll(elementData, j)
    {
        // Read tower element data
        scalar axialDistance = elementData[j][0];
        scalar height = elementData[j][1];
        scalar diameter = elementData[j][2];
        
        // Set sizes for actuatorLineSource elementGeometry lists
        elementGeometry[j].setSize(6);
        elementGeometry[j][0].setSize(3);
        elementGeometry[j][1].setSize(3);
        elementGeometry[j][2].setSize(1);
        elementGeometry[j][3].setSize(3);
        elementGeometry[j][4].setSize(1);
        elementGeometry[j][5].setSize(1);
        
        // Create geometry point for AL source at origin
        vector point = origin_;
        // Move along turbine axis
        point += axialDistance*axis_;
        // Move along tower axis according to height
        point += height*towerAxis;
        
        elementGeometry[j][0][0] = point.x(); // x location of geom point
        elementGeometry[j][0][1] = point.y(); // y location of geom point
        elementGeometry[j][0][2] = point.z(); // z location of geom point
        
        // Set span directions
        elementGeometry[j][1][0] = towerAxis.x(); // x component of span
        elementGeometry[j][1][1] = towerAxis.y(); // y component of span
        elementGeometry[j][1][2] = towerAxis.z(); // z component of span
        
        // Set chord length
        elementGeometry[j][2][0] = diameter;
        
        // Set chord reference direction
        elementGeometry[j][3][0] = freeStreamDirection_.x();
        elementGeometry[j][3][1] = freeStreamDirection_.y();
        elementGeometry[j][3][2] = freeStreamDirection_.z();
        
        // Set chord mount
        elementGeometry[j][4][0] = 0.25;
        
        // Set pitch
        elementGeometry[j][5][0] = 0.0;
    }
    
    towerSubDict.add("elementGeometry", elementGeometry);
    towerSubDict.add("initialVelocities", initialVelocities);
    towerSubDict.add("fieldNames", coeffs_.lookup("fieldNames"));
    towerSubDict.add("profileData", profileData_);
    towerSubDict.add("tipEffect", tipEffect_);
    towerSubDict.add("freeStreamVelocity", freeStreamVelocity_);
        
    dictionary dict;
    dict.add("actuatorLineSourceCoeffs", towerSubDict);
    dict.add("type", "actuatorLineSource");
    dict.add("active", dict_.lookup("active"));
    dict.add("selectionMode", dict_.lookup("selectionMode"));
    dict.add("cellSet", dict_.lookup("cellSet"));
    
    actuatorLineSource* tower = new actuatorLineSource
    (
        "tower", 
        "actuatorLineSource", 
        dict, 
        mesh_
    );
    
    tower_.set(tower);
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::fv::axialFlowTurbineALSource::axialFlowTurbineALSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    turbineALSource(name, modelType, dict, mesh),
    hasHub_(false),
    hasTower_(false),
    hasNacelle_(false),
    verticalDirection_
    (
        coeffs_.lookupOrDefault("verticalDirection", vector(0, 0, 1))
    )
{
    read(dict);
    createCoordinateSystem();
    createBlades();
    if (hasHub_) createHub();
    if (hasTower_) createTower();
    if (hasNacelle_) createNacelle();
    createOutputFile();
    
    if (debug)
    {
        Info<< "axialFlowTurbineALSource created at time = " << time_.value()
            << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::axialFlowTurbineALSource::~axialFlowTurbineALSource()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::axialFlowTurbineALSource::rotate()
{
    // Update tip speed ratio and omega
    scalar t = time_.value();
    tipSpeedRatio_ = meanTSR_ + tsrAmplitude_
                   *cos(nBlades_/rotorRadius_*(meanTSR_*t - tsrPhase_));
    omega_ = tipSpeedRatio_*mag(freeStreamVelocity_)/rotorRadius_;

    scalar deltaT = time_.deltaT().value();
    scalar radians = omega_*deltaT;
    
    forAll(blades_, i)
    {
        blades_[i].rotate(origin_, axis_, radians);
        blades_[i].setSpeed(origin_, axis_, omega_);
    }
    
    if (hasHub_)
    {
        hub_->rotate(origin_, axis_, radians);
        hub_->setSpeed(origin_, axis_, omega_);
    }
    
    if (debug)
    {
        Info<< "Rotating " << name_ << " " << radians << " radians" 
            << endl << endl;
    }
    angleDeg_ += radians*180.0/Foam::constant::mathematical::pi;
    lastRotationTime_ = time_.value();
}


void Foam::fv::axialFlowTurbineALSource::addSup
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
    
    // Create local moment vector
    vector moment(vector::zero);

    // Add source for blade actuator lines
    forAll(blades_, i)
    {
        blades_[i].addSup(eqn, fieldI);
        forceField_ += blades_[i].forceField();
        force_ += blades_[i].force();
        moment += blades_[i].moment(origin_);
    }
    
    if (hasHub_)
    {
        // Add source for hub actuator line
        hub_->addSup(eqn, fieldI);
        forceField_ += hub_->forceField();
        force_ += hub_->force();
        moment += hub_->moment(origin_);
    }
    
    if (hasTower_)
    {
        // Add source for tower actuator line
        tower_->addSup(eqn, fieldI);
        forceField_ += tower_->forceField();
        if (includeTowerDrag_) force_ += tower_->force();
    }
    
    if (hasNacelle_)
    {
        // Add source for tower actuator line
        nacelle_->addSup(eqn, fieldI);
        forceField_ += nacelle_->forceField();
        if (includeNacelleDrag_) force_ += nacelle_->force();
    }
    
    // Torque is the projection of the moment from all blades on the axis
    torque_ = moment & axis_;
    Info<< "Azimuthal angle (degrees) of " << name_ << ": " << angleDeg_ 
        << endl;
    Info<< "Torque (per unit density) from " << name_ << ": " << torque_ 
        << endl;
    
    torqueCoefficient_ = torque_/(0.5*frontalArea_*rotorRadius_
                       * magSqr(freeStreamVelocity_));
    powerCoefficient_ = torqueCoefficient_*tipSpeedRatio_;
    dragCoefficient_ = force_ & freeStreamDirection_
                     / (0.5*frontalArea_*magSqr(freeStreamVelocity_));
                             
    Info<< "Power coefficient from " << name_ << ": " << powerCoefficient_
        << endl << endl;
        
    // Write performance data -- note this will write multiples if there are
    // multiple PIMPLE loops
    if (Pstream::master()) writePerf();
}


void Foam::fv::axialFlowTurbineALSource::addSup
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
    
    // Create local moment vector
    vector moment(vector::zero);

    // Add source for blade actuator lines
    forAll(blades_, i)
    {
        blades_[i].addSup(rho, eqn, fieldI);
        forceField_ += blades_[i].forceField();
        force_ += blades_[i].force();
        moment += blades_[i].moment(origin_);
    }
    
    if (hasHub_)
    {
        // Add source for hub actuator line
        hub_->addSup(rho, eqn, fieldI);
        forceField_ += hub_->forceField();
        force_ += hub_->force();
        moment += hub_->moment(origin_);
    }
    
    if (hasTower_)
    {
        // Add source for tower actuator line
        tower_->addSup(rho, eqn, fieldI);
        forceField_ += tower_->forceField();
        if (includeTowerDrag_) force_ += tower_->force();
    }
    
    if (hasNacelle_)
    {
        // Add source for tower actuator line
        nacelle_->addSup(rho, eqn, fieldI);
        forceField_ += nacelle_->forceField();
        if (includeNacelleDrag_) force_ += nacelle_->force();
    }
    
    // Torque is the projection of the moment from all blades on the axis
    torque_ = moment & axis_;
    Info<< "Azimuthal angle (degrees) of " << name_ << ": " << angleDeg_ 
        << endl;
    Info<< "Torque from " << name_ << ": " << torque_ 
        << endl;
    
    scalar rhoRef;
    coeffs_.lookup("rhoRef") >> rhoRef;
    torqueCoefficient_ = torque_/(0.5*rhoRef*frontalArea_*rotorRadius_
                       * magSqr(freeStreamVelocity_));
    powerCoefficient_ = torqueCoefficient_*tipSpeedRatio_;
    dragCoefficient_ = force_ & freeStreamDirection_
                     / (0.5*rhoRef*frontalArea_*magSqr(freeStreamVelocity_));
                             
    Info<< "Power coefficient from " << name_ << ": " << powerCoefficient_
        << endl << endl;
        
    // Write performance data -- note this will write multiples if there are
    // multiple PIMPLE loops
    if (Pstream::master()) writePerf();
}


void Foam::fv::axialFlowTurbineALSource::addSup
(
    fvMatrix<scalar>& eqn,
    const label fieldI
)
{
    // Rotate the turbine if time value has changed
    if (time_.value() != lastRotationTime_)
    {
        rotate();
    }

    // Add scalar source term from blades
    forAll(blades_, i)
    {
        blades_[i].addSup(eqn, fieldI);
    }
    
    if (hasHub_)
    {
        // Add source for hub actuator line
        hub_->addSup(eqn, fieldI);
    }
    
    if (hasTower_)
    {
        // Add source for tower actuator line
        tower_->addSup(eqn, fieldI);
    }
    
    if (hasNacelle_)
    {
        // Add source for nacelle actuator line
        nacelle_->addSup(eqn, fieldI);
    }
}


void Foam::fv::axialFlowTurbineALSource::printCoeffs() const
{
    Info<< "Number of blades: " << nBlades_ << endl;
}


bool Foam::fv::axialFlowTurbineALSource::read(const dictionary& dict)
{
    if (option::read(dict))
    {
        turbineALSource::read(dict);
        
        // Get hub information
        hubDict_ = coeffs_.subOrEmptyDict("hub");
        if (hubDict_.keys().size() > 0) hasHub_ = true;
        
        // Get tower information
        towerDict_ = coeffs_.subOrEmptyDict("tower");
        if (towerDict_.keys().size() > 0) hasTower_ = true;
        includeTowerDrag_ = towerDict_.lookupOrDefault
        (
            "includeInTotalDrag",
            false
        );
        
        // Get nacelle information
        nacelleDict_ = coeffs_.subOrEmptyDict("nacelle");
        if (nacelleDict_.keys().size() > 0) hasNacelle_ = true;
        includeNacelleDrag_ = nacelleDict_.lookupOrDefault
        (
            "includeInTotalDrag",
            false
        );
        
        if (debug)
        {
            Info<< "Debugging on" << endl;
            Info<< "Axial-flow turbine properties:" << endl;
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
