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

#include "crossFlowTurbineALSource.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMatrices.H"
#include "geometricOneField.H"
#include "syncTools.H"
#include "unitConversion.H"

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

void Foam::fv::crossFlowTurbineALSource::createCoordinateSystem()
{
    // Construct the local rotor coordinate system
    freeStreamDirection_ = freeStreamVelocity_/mag(freeStreamVelocity_);
    radialDirection_ = axis_^freeStreamDirection_;
    radialDirection_ = radialDirection_/mag(radialDirection_);
    // Make sure axis is a unit vector
    axis_ /= mag(axis_);
}


void Foam::fv::crossFlowTurbineALSource::createBlades()
{
    int nBlades = nBlades_;
    blades_.setSize(nBlades);
    int nElements;
    List<List<scalar> > elementData;
    List<List<scalar> > profileData;
    word modelType = "actuatorLineSource";
    List<scalar> frontalAreas(nBlades); // Frontal area from each blade

    forAll(blades_, i)
    {
        word bladeName = bladeNames_[i];
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
        forAll(elementData, j)
        {
            // Read CFTAL dict data
            scalar axialDistance = elementData[j][0];
            scalar radius = elementData[j][1];
            scalar azimuthDegrees = elementData[j][2] + azimuthalOffset;
            scalar azimuthRadians = degToRad(azimuthDegrees);
            scalar chordLength = elementData[j][3];
            scalar chordMount = elementData[j][4];
            scalar pitch = elementData[j][5];

            // Compute frontal area contribution from this geometry segment
            if (j > 0)
            {
                scalar deltaAxial = axialDistance - elementData[j-1][0];
                scalar meanRadius = (radius + elementData[j-1][1])/2;
                frontalArea += mag(deltaAxial*meanRadius);
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
            // Move along axis
            point += axialDistance*axis_;
            // Move along chord according to chordMount
            scalar chordDisplacement = (chordMount - 0.25)*chordLength;
            point -= chordDisplacement*freeStreamDirection_;
            // Move along radial direction
            point += radius*radialDirection_;
            // Set initial velocity of quarter chord
            scalar radiusCorr = sqrt(magSqr((chordMount - 0.25)*chordLength)
                                     + magSqr(radius));
            vector initialVelocity = -freeStreamDirection_*omega_*radiusCorr;
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
            elementGeometry[j][1][0] = axis_.x(); // x component of span dir
            elementGeometry[j][1][1] = axis_.y(); // y component of span dir
            elementGeometry[j][1][2] = axis_.z(); // z component of span dir

            // Set chord length
            elementGeometry[j][2][0] = chordLength;

            // Set chord reference direction
            vector chordDirection = -freeStreamDirection_;
            rotateVector(chordDirection, vector::zero, axis_, azimuthRadians);
            elementGeometry[j][3][0] = chordDirection.x();
            elementGeometry[j][3][1] = chordDirection.y();
            elementGeometry[j][3][2] = chordDirection.z();

            // Set chord mount
            elementGeometry[j][4][0] = chordMount;

            // Set pitch
            elementGeometry[j][5][0] = pitch;
        }

        // Add frontal area to list
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
        bladeSubDict.add
        (
            "addedMass",
            coeffs_.lookupOrDefault("addedMass", false)
        );
        bladeSubDict.add
        (
            "velocitySampleRadius",
            coeffs_.lookupOrDefault("velocitySampleRadius", 0.0)
        );
        bladeSubDict.add
        (
            "nVelocitySamples",
            coeffs_.lookupOrDefault("nVelocitySamples", 20)
        );
        bladeSubDict.add("selectionMode", coeffs_.lookup("selectionMode"));
        bladeSubDict.add("cellSet", coeffs_.lookup("cellSet"));

        // Lookup or create flowCurvature subDict
        dictionary fcDict = coeffs_.subOrEmptyDict("flowCurvature");
        fcDict.lookupOrAddDefault("active", true);
        word defaultFCModel = "Goude";
        fcDict.lookupOrAddDefault
        (
            "flowCurvatureModel",
            defaultFCModel
        );
        bladeSubDict.add("flowCurvature", fcDict);

        // Do not write force from individual actuator line unless specified
        bladeSubDict.lookupOrAddDefault("writeForceField", false);

        dictionary dict;
        dict.add("actuatorLineSourceCoeffs", bladeSubDict);
        dict.add("type", "actuatorLineSource");
        dict.add("active", dict_.lookup("active"));

        actuatorLineSource* blade = new actuatorLineSource
        (
            name_ + "." + bladeName,
            modelType,
            dict,
            mesh_
        );

        blades_.set(i, blade);
    }

    // Frontal area is twice the maximum blade frontal area
    frontalArea_ = 2*max(frontalAreas);
    Info<< "Frontal area of " << name_ << ": " << frontalArea_ << endl;
}


void Foam::fv::crossFlowTurbineALSource::createStruts()
{
    int nStruts = strutsDict_.keys().size();
    struts_.setSize(nStruts);
    int nElements;
    List<List<scalar> > elementData;
    word modelType = "actuatorLineSource";
    List<word> strutNames = strutsDict_.toc();

    forAll(struts_, i)
    {
        word strutName = strutNames[i];
        // Create dictionary items for this strut
        dictionary strutSubDict = strutsDict_.subDict(strutName);
        strutSubDict.lookup("nElements") >> nElements;
        strutSubDict.lookup("elementData") >> elementData;
        scalar azimuthalOffset = strutSubDict.lookupOrDefault
        (
            "azimuthalOffset",
            0.0
        );

        strutSubDict.add("freeStreamVelocity", freeStreamVelocity_);
        strutSubDict.add("fieldNames", coeffs_.lookup("fieldNames"));
        strutSubDict.add("profileData", profileData_);

        if (debug)
        {
            Info<< "Creating actuator line strut " << strutName << endl;
            Info<< "Strut has " << nElements << " elements" << endl;
            Info<< "Element data:" << endl;
            Info<< elementData << endl << endl;
        }

        // Convert element data into actuator line element geometry
        label nGeomPoints = elementData.size();
        List<List<List<scalar> > > elementGeometry(nGeomPoints);
        List<vector> initialVelocities(nGeomPoints, vector::one);

        forAll(elementData, j)
        {
            // Read CFTAL dict data
            scalar axialDistance = elementData[j][0];
            scalar radius = elementData[j][1];
            scalar azimuthDegrees = elementData[j][2] + azimuthalOffset;
            scalar azimuthRadians = degToRad(azimuthDegrees);
            scalar chordLength = elementData[j][3];
            scalar chordMount = elementData[j][4];
            scalar pitch = elementData[j][5];

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
            // Move along chord according to chordMount
            scalar chordDisplacement = (chordMount - 0.25)*chordLength;
            point -= chordDisplacement*freeStreamDirection_;
            // Move along radial direction
            point += radius*radialDirection_;
            // Set initial velocity of quarter chord
            scalar radiusCorr = sqrt(magSqr((chordMount - 0.25)*chordLength)
                                     + magSqr(radius));
            vector initialVelocity = -freeStreamDirection_*omega_*radiusCorr;
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

            // Set span directions for AL source (in radial direction)
            vector spanDirection = radialDirection_;
            rotateVector(spanDirection, vector::zero, axis_, azimuthRadians);
            elementGeometry[j][1][0] = spanDirection.x();
            elementGeometry[j][1][1] = spanDirection.y();
            elementGeometry[j][1][2] = spanDirection.z();

            // Set chord length
            elementGeometry[j][2][0] = chordLength;

            // Set chord reference direction
            vector chordDirection = -freeStreamDirection_;
            rotateVector(chordDirection, vector::zero, axis_, azimuthRadians);
            elementGeometry[j][3][0] = chordDirection.x();
            elementGeometry[j][3][1] = chordDirection.y();
            elementGeometry[j][3][2] = chordDirection.z();

            // Set chord mount
            elementGeometry[j][4][0] = chordMount;

            // Set pitch
            elementGeometry[j][5][0] = pitch;
        }

        if (debug)
        {
            Info<< "Converted element geometry:" << endl << elementGeometry
                << endl;
        }

        strutSubDict.add("elementGeometry", elementGeometry);
        strutSubDict.add("initialVelocities", initialVelocities);
        strutSubDict.add("selectionMode", coeffs_.lookup("selectionMode"));
        strutSubDict.add("cellSet", coeffs_.lookup("cellSet"));

        // Do not write force from individual actuator line unless specified
        strutSubDict.lookupOrAddDefault("writeForceField", false);

        dictionary dict;
        dict.add("actuatorLineSourceCoeffs", strutSubDict);
        dict.add("type", "actuatorLineSource");
        dict.add("active", dict_.lookup("active"));

        actuatorLineSource* strut = new actuatorLineSource
        (
            name_ + "." + strutName,
            modelType,
            dict,
            mesh_
        );

        struts_.set(i, strut);
    }
}


void Foam::fv::crossFlowTurbineALSource::createShaft()
{
    List<List<scalar> > elementData;
    dictionary shaftSubDict = shaftDict_;

    shaftDict_.lookup("elementData") >> elementData;

    // Convert element data into actuator line element geometry
    label nGeomPoints = elementData.size();
    List<List<List<scalar> > > elementGeometry(nGeomPoints);
    List<vector> initialVelocities(nGeomPoints, vector::zero);

    forAll(elementData, j)
    {
        // Read shaft element data
        scalar axialDistance = elementData[j][0];
        scalar diameter = elementData[j][1];

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

        elementGeometry[j][0][0] = point.x(); // x location of geom point
        elementGeometry[j][0][1] = point.y(); // y location of geom point
        elementGeometry[j][0][2] = point.z(); // z location of geom point

        // Set span directions
        elementGeometry[j][1][0] = axis_.x(); // x component of span direction
        elementGeometry[j][1][1] = axis_.y(); // y component of span direction
        elementGeometry[j][1][2] = axis_.z(); // z component of span direction

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

    shaftSubDict.add("elementGeometry", elementGeometry);
    shaftSubDict.add("initialVelocities", initialVelocities);
    shaftSubDict.add("fieldNames", coeffs_.lookup("fieldNames"));
    shaftSubDict.add("profileData", profileData_);
    shaftSubDict.add("freeStreamVelocity", freeStreamVelocity_);
    shaftSubDict.add("selectionMode", coeffs_.lookup("selectionMode"));
    shaftSubDict.add("cellSet", coeffs_.lookup("cellSet"));

    // Do not write force from individual actuator line unless specified
    shaftSubDict.lookupOrAddDefault("writeForceField", false);

    dictionary dict;
    dict.add("actuatorLineSourceCoeffs", shaftSubDict);
    dict.add("type", "actuatorLineSource");
    dict.add("active", dict_.lookup("active"));

    actuatorLineSource* shaft = new actuatorLineSource
    (
        name_ + ".shaft",
        "actuatorLineSource",
        dict,
        mesh_
    );

    shaft_.set(shaft);
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
    turbineALSource(name, modelType, dict, mesh),
    hasStruts_(false),
    hasShaft_(false)
{
    read(dict);
    createCoordinateSystem();
    createBlades();
    if (hasStruts_)
    {
        createStruts();
    }
    if (hasShaft_)
    {
        createShaft();
    }
    createOutputFile();

    // Rotate turbine to azimuthalOffset if necessary
    scalar azimuthalOffset = coeffs_.lookupOrDefault("azimuthalOffset", 0.0);
    rotate(degToRad(azimuthalOffset));

    if (debug)
    {
        Info<< "crossFlowTurbineALSource created at time = " << time_.value()
            << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::crossFlowTurbineALSource::~crossFlowTurbineALSource()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::crossFlowTurbineALSource::rotate(scalar radians)
{
    if (debug)
    {
        Info<< "Rotating " << name_ << " " << radians << " radians"
            << endl << endl;
    }

    forAll(blades_, i)
    {
        blades_[i].rotate(origin_, axis_, radians);
        blades_[i].setSpeed(origin_, axis_, omega_);
    }

    if (hasStruts_)
    {
        forAll(struts_, i)
        {
            struts_[i].rotate(origin_, axis_, radians);
            struts_[i].setSpeed(origin_, axis_, omega_);
        }
    }

    if (hasShaft_)
    {
        shaft_->rotate(origin_, axis_, radians);
        shaft_->setSpeed(origin_, axis_, omega_);
    }
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
    forceField_ *= dimensionedScalar("zero", forceField_.dimensions(), 0.0);
    force_ *= 0;

    // Check dimensions of force field and correct if necessary
    if (forceField_.dimensions() != eqn.dimensions()/dimVolume)
    {
        forceField_.dimensions().reset(eqn.dimensions()/dimVolume);
    }

    // Create local moment vector
    vector moment(vector::zero);

    // Add source for blade actuator lines
    forAll(blades_, i)
    {
        blades_[i].addSup(eqn, fieldI);
        forceField_ += blades_[i].forceField();
        Info<< "Added blade" << endl;
        force_ += blades_[i].force();
        bladeMoments_[i] = blades_[i].moment(origin_);
        moment += bladeMoments_[i];
    }

    if (hasStruts_)
    {
        // Add source for strut actuator lines
        forAll(struts_, i)
        {
            struts_[i].addSup(eqn, fieldI);
            forceField_ += struts_[i].forceField();
            force_ += struts_[i].force();
            moment += struts_[i].moment(origin_);
        }
    }

    if (hasShaft_)
    {
        // Add source for shaft actuator line
        shaft_->addSup(eqn, fieldI);
        forceField_ += shaft_->forceField();
        force_ += shaft_->force();
        moment += shaft_->moment(origin_);
    }

    // Torque is the projection of the moment from all blades on the axis
    torque_ = moment & axis_;

    torqueCoefficient_ = torque_/(0.5*frontalArea_*rotorRadius_
                       * magSqr(freeStreamVelocity_));
    powerCoefficient_ = torqueCoefficient_*tipSpeedRatio_;
    dragCoefficient_ = force_ & freeStreamDirection_
                     / (0.5*frontalArea_*magSqr(freeStreamVelocity_));

    // Print performance to terminal
    printPerf();

    // Write performance data -- note this will write multiples if there are
    // multiple PIMPLE loops
    if (Pstream::master())
    {
        writePerf();
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
    forceField_ *= dimensionedScalar("zero", forceField_.dimensions(), 0.0);
    force_ *= 0;

    // Create local moment vector
    vector moment(vector::zero);

    // Check dimensions of force field and correct if necessary
    if (forceField_.dimensions() != eqn.dimensions()/dimVolume)
    {
        forceField_.dimensions().reset(eqn.dimensions()/dimVolume);
    }

    // Add source for blade actuator lines
    forAll(blades_, i)
    {
        blades_[i].addSup(rho, eqn, fieldI);
        forceField_ += blades_[i].forceField();
        force_ += blades_[i].force();
        bladeMoments_[i] = blades_[i].moment(origin_);
        moment += bladeMoments_[i];
    }

    if (hasStruts_)
    {
        // Add source for strut actuator lines
        forAll(struts_, i)
        {
            struts_[i].addSup(rho, eqn, fieldI);
            forceField_ += struts_[i].forceField();
            force_ += struts_[i].force();
            moment += struts_[i].moment(origin_);
        }
    }

    if (hasShaft_)
    {
        // Add source for shaft actuator line
        shaft_->addSup(rho, eqn, fieldI);
        forceField_ += shaft_->forceField();
        force_ += shaft_->force();
        moment += shaft_->moment(origin_);
    }

    // Torque is the projection of the moment from all blades on the axis
    torque_ = moment & axis_;

    scalar rhoRef;
    coeffs_.lookup("rhoRef") >> rhoRef;
    torqueCoefficient_ = torque_/(0.5*rhoRef*frontalArea_*rotorRadius_
                       * magSqr(freeStreamVelocity_));
    powerCoefficient_ = torqueCoefficient_*tipSpeedRatio_;
    dragCoefficient_ = force_ & freeStreamDirection_
                     / (0.5*rhoRef*frontalArea_*magSqr(freeStreamVelocity_));

    // Print performance to terminal
    printPerf();

    // Write performance data -- note this will write multiples if there are
    // multiple PIMPLE loops
    if (Pstream::master())
    {
        writePerf();
    }
}


void Foam::fv::crossFlowTurbineALSource::addSup
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

    if (hasStruts_)
    {
        // Add source for strut actuator lines
        forAll(struts_, i)
        {
            struts_[i].addSup(eqn, fieldI);
        }
    }

    if (hasShaft_)
    {
        // Add source for shaft actuator line
        shaft_->addSup(eqn, fieldI);
    }
}


void Foam::fv::crossFlowTurbineALSource::printCoeffs() const
{
    Info<< "Number of blades: " << nBlades_ << endl;
}


bool Foam::fv::crossFlowTurbineALSource::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    {
        turbineALSource::read(dict);

        // Get struts information
        strutsDict_ = coeffs_.subOrEmptyDict("struts");
        if (strutsDict_.keys().size() > 0)
        {
            hasStruts_ = true;
        }

        // Get shaft information
        shaftDict_ = coeffs_.subOrEmptyDict("shaft");
        if (shaftDict_.keys().size() > 0)
        {
            hasShaft_ = true;
        }

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
