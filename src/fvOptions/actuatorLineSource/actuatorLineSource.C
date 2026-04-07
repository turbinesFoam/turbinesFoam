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

#include "actuatorLineSource.H"
#include "unitConversion.H"
#include "addToRunTimeSelectionTable.H"
#include "vector.H"
#include "fvMatrices.H"
#include "geometricOneField.H"
#include "syncTools.H"
#include "simpleMatrix.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(actuatorLineSource, 0);
    addToRunTimeSelectionTable
    (
        option,
        actuatorLineSource,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

bool Foam::fv::actuatorLineSource::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    {

        coeffs_.lookup("fieldNames") >> fieldNames_;
        applied_.setSize(fieldNames_.size(), false);

        // Look up information in dictionary
        coeffs_.lookup("elementProfiles") >> elementProfiles_;
        profileData_ = coeffs_.subDict("profileData");
        coeffs_.lookup("elementGeometry") >> elementGeometry_;
        coeffs_.lookup("nElements") >> nElements_;
        coeffs_.lookup("freeStreamVelocity") >> freeStreamVelocity_;
        freeStreamDirection_ = freeStreamVelocity_/mag(freeStreamVelocity_);
        endEffectsActive_ = coeffs_.lookupOrDefault("endEffects", false);

        // Read harmonic pitching parameters if present
        dictionary pitchDict = coeffs_.subOrEmptyDict("harmonicPitching");
        harmonicPitchingActive_ = pitchDict.lookupOrDefault("active", false);
        reducedFreq_ = pitchDict.lookupOrDefault("reducedFreq", 0.0);
        pitchAmplitude_ = pitchDict.lookupOrDefault("amplitude", 0.0);

        // Read option for writing forceField
        bool writeForceField = coeffs_.lookupOrDefault
        (
            "writeForceField",
            true
        );
        if (not writeForceField)
        {
            forceField_.writeOpt() = IOobject::NO_WRITE;
        }

        if (debug)
        {
            Info<< "Debugging for actuatorLineSource on" << endl;
            printCoeffs();
        }

        return true;
    }
    else
    {
        return false;
    }
}


void Foam::fv::actuatorLineSource::createOutputFile()
{
    fileName dir;

    if (Pstream::parRun())
    {
        dir = mesh_.time().path()/"../postProcessing/actuatorLines"
            / mesh_.time().timeName();
    }
    else
    {
        dir = mesh_.time().path()/"postProcessing/actuatorLines"
            / mesh_.time().timeName();
    }

    if (not isDir(dir))
    {
        mkDir(dir);
    }

    outputFile_ = new OFstream(dir/name_ + ".csv");

    *outputFile_<< "time,x,y,z,rel_vel_mag,alpha_deg,alpha_geom_deg,cl,cd,cm"
                << endl;
}


void Foam::fv::actuatorLineSource::createElements()
{
    elements_.setSize(nElements_);

    label nGeometryPoints = elementGeometry_.size();
    //label nGeometrySegments = nGeometryPoints - 1;
    label nGeometrySegments = nGeometryPoints;
    label nElementsPerSegment = nElements_/nGeometrySegments;
    if (nElements_ % nGeometrySegments)
    {
        // Need to have integer number of elements per geometry segment
        FatalErrorIn("void actuatorLineSource::createElements()")
            << "Number of actuator line elements must be multiple of the "
            << "number of actuator line geometry segments"
            << abort(FatalError);
    }
    List<vector> points(nGeometryPoints);
    List<vector> spanDirs(nGeometryPoints);
    List<scalar> chordLengths(nGeometryPoints);
    List<scalar> spanLengths(nGeometrySegments);
    List<vector> chordRefDirs(nGeometryPoints);
    List<scalar> pitches(nGeometryPoints);
    List<scalar> chordMounts(nGeometryPoints);
    totalLength_ = 0.0;
    chordLength_ = 0.0;

    forAll(points, i)
    {
        // Extract geometry point
        scalar x = elementGeometry_[i][0][0];
        scalar y = elementGeometry_[i][0][1];
        scalar z = elementGeometry_[i][0][2];
        points[i] = vector(x, y, z);
    }

    // Modification in the way the individual BEM elements are calculated, to ensure that the whole blade
    // radius is accounted for
    // Store blade root and tip locations for distance calculations
    vector rootLocation;
    rootLocation[0] = elementGeometry_[0][6][0]; // Initialization of rootLocation at rotor center
    rootLocation[1] = elementGeometry_[0][6][1];
    rootLocation[2] = elementGeometry_[0][6][2];
    vector tipLocation;
    tipLocation[0] = elementGeometry_[nGeometryPoints-1][7][0]; // Initialization of tipLocation at rotor tip
    tipLocation[1] = elementGeometry_[nGeometryPoints-1][7][1];
    tipLocation[2] = elementGeometry_[nGeometryPoints-1][7][2];
    
    List<vector> BEMpoint1(nGeometryPoints); // First point of BEM element
    List<vector> BEMpoint2(nGeometryPoints); // Second point of BEM element
    BEMpoint1[0] = rootLocation; // Initialization of BEMpoint1 for first BEM element at rotor center
 
    forAll(points, i)
    {
        if (i != nGeometryPoints - 1)
        {
            // Linearly interpolate position of current BEM endpoint
            BEMpoint2[i] = (points[i] + points[i + 1])/2;
            spanLengths[i] = mag(BEMpoint2[i] - BEMpoint1[i]);
            totalLength_ += spanLengths[i];
            BEMpoint1[i+1] = BEMpoint2[i]; // Update BEMpoint1 for next iteration
        } else {
            // Last BEm endpoint equal to blade tip point     
            BEMpoint2[i] = tipLocation; // Initialization of BEMpoint2 at rotor tip
            spanLengths[i] = mag(BEMpoint2[i] - BEMpoint1[i]);
            totalLength_ += spanLengths[i];
        }
        // Read span direction
        scalar x = elementGeometry_[i][1][0];
        scalar y = elementGeometry_[i][1][1];
        scalar z = elementGeometry_[i][1][2];
        spanDirs[i] = vector(x, y, z);
        // Read chord length
        chordLengths[i] = elementGeometry_[i][2][0];
        chordLength_ += chordLengths[i];
        // Read chord ref dir
        x = elementGeometry_[i][3][0];
        y = elementGeometry_[i][3][1];
        z = elementGeometry_[i][3][2];
        chordRefDirs[i] = vector(x, y, z);
        // Read chord mount
        chordMounts[i] = elementGeometry_[i][4][0];
        // Read pitch
        pitches[i] = elementGeometry_[i][5][0];
    }

    // Compute average chord length
    chordLength_ /= nGeometryPoints;

    // Compute aspect ratio
    aspectRatio_ = totalLength_/chordLength_;

    // Lookup initial element velocities if present
    List<vector> initialVelocities(nGeometryPoints, vector::zero);
    coeffs_.readIfPresent("initialVelocities", initialVelocities);

    if (debug)
    {
        Info<< "Total length: " << totalLength_ << endl;
        Info<< "Elements per geometry segment: " << nElementsPerSegment
            << endl;
        Info<< "Points:" << endl << points << endl;
        Info<< "Span directions:" << endl << spanDirs << endl;
        Info<< "Span lengths: " << endl << spanLengths << endl;
        Info<< "Chord lengths:" << endl << chordLengths << endl;
        Info<< "Pitches:" << endl << pitches << endl;
        Info<< "Root location: " << rootLocation << endl;
        Info<< "Tip location: " << tipLocation << endl;
    }

    forAll(elements_, i)
    {
        std::stringstream ss;
        ss << i;
        string str = ss.str();
        const word name = name_ + ".element" + str;

        // Actuator point geometry to be calculated from elementGeometry
        label geometrySegmentIndex = i/nElementsPerSegment;
        label pointIndex = i % nElementsPerSegment;
        label elementProfileIndex = i*elementProfiles_.size()/nElements_;
        word profileName = elementProfiles_[elementProfileIndex];
        vector position;
        scalar chordLength;
        vector chordDirection;
        vector chordRefDirection;
        scalar spanLength = spanLengths[geometrySegmentIndex];
        spanLength /= nElementsPerSegment;
        vector spanDirection;
        scalar pitch;
        scalar chordMount;
        vector initialVelocity;

        // Linearly interpolate position
        vector point1 = BEMpoint1[geometrySegmentIndex];
        vector point2 = BEMpoint2[geometrySegmentIndex];
        vector segment = point2 - point1;
        position = point1
                 + segment/nElementsPerSegment*pointIndex
                 + segment/nElementsPerSegment/2;
				 
        // The geometry characteristics are currently constant on each blade element
        // Read geometric chordLength
        chordLength = chordLengths[geometrySegmentIndex];

        // Read geometric spanDirection
        spanDirection = spanDirs[geometrySegmentIndex];

        // Read geometric section pitch
        pitch = pitches[geometrySegmentIndex];

        // Read geometric chord mount
        chordMount = chordMounts[geometrySegmentIndex];

        // Read geometric element velocity
        initialVelocity = initialVelocities[geometrySegmentIndex];

        // Read geometric chordDirection
        chordDirection = chordRefDirs[geometrySegmentIndex];
        // Chord reference direction (before pitching)
        chordRefDirection = chordDirection;
        
        // Calculate nondimensional root distance
        scalar rootDistance = mag(position - rootLocation)/totalLength_;

        // Create a dictionary for this actuatorLineElement
        dictionary dict;
        dict.add("position", position);
        dictionary profileDataDict = profileData_.subDict(profileName);
        dict.add("profileData", profileDataDict);
        dict.add("profileName", profileName);
        dict.add("pitch", pitch);
        dict.add("chordLength", chordLength);
        dict.add("chordDirection", chordDirection);
        dict.add("chordRefDirection", chordRefDirection);
        dict.add("spanLength", spanLength);
        dict.add("spanDirection", spanDirection);
        dict.add("freeStreamVelocity", freeStreamVelocity_);
        dict.add("chordMount", chordMount);
        dict.add("rootDistance", rootDistance);
        dict.add("addedMass", coeffs_.lookupOrDefault("addedMass", false));
        dict.add
        (
            "velocitySampleRadius",
            coeffs_.lookupOrDefault("velocitySampleRadius", 0.0)
        );
        dict.add
        (
            "nVelocitySamples",
            coeffs_.lookupOrDefault("nVelocitySamples", 20)
        );
        if (coeffs_.found("dynamicStall"))
        {
            dictionary dsDict = coeffs_.subDict("dynamicStall");
            dsDict.add("chordLength", chordLength);
            dict.add("dynamicStall", dsDict);
        }
        dictionary fcDict = coeffs_.subOrEmptyDict("flowCurvature");
        dict.add("flowCurvature", fcDict);
        bool writeElementPerf
        (
            coeffs_.lookupOrDefault("writeElementPerf", false)
        );
        dict.add("writePerf", writeElementPerf);

        if (debug)
        {
            Info<< "Creating actuatorLineElement: " << name << endl;
            Info<< "Geometry segment index: " << geometrySegmentIndex << endl;
            Info<< "Position: " << position << endl;
            Info<< "Chord length: " << chordLength << endl;
            Info<< "Chord direction (before pitching): " << chordDirection
                << endl;
            Info<< "Pitch (degrees): " << pitch << endl;
            Info<< "Span length: " << spanLength << endl;
            Info<< "Span direction: " << spanDirection << endl;
            Info<< "Profile name index: " << elementProfileIndex << endl;
            Info<< "Profile name: " << profileName << endl;
            Info<< "writePerf: " << writeElementPerf << endl;
            Info<< "Root distance (nondimensional): " << rootDistance << endl;
        }

        actuatorLineElement* element = new actuatorLineElement
        (
            name, dict, mesh_
        );
        elements_.set(i, element);
        pitch = Foam::degToRad(pitch);
        elements_[i].pitch(pitch);
        elements_[i].setVelocity(initialVelocity);
    }
}


void Foam::fv::actuatorLineSource::writePerf()
{
    scalar time = mesh_.time().value();
    scalar totalArea = 0.0;
    scalar x = 0.0;
    scalar y = 0.0;
    scalar z = 0.0;
    scalar relVelMag = 0.0;
    scalar alphaDeg = 0.0;
    scalar alphaGeom = 0.0;
    scalar cl = 0.0;
    scalar cd = 0.0;
    scalar cm = 0.0;

    forAll(elements_, i)
    {
        scalar area = elements_[i].chordLength()*elements_[i].spanLength();
        totalArea += area;
        vector pos = elements_[i].position();
        x += pos[0]; y += pos[1]; z += pos[2];
        relVelMag += mag(elements_[i].relativeVelocity())*area;
        alphaDeg += elements_[i].angleOfAttack()*area;
        alphaGeom += elements_[i].angleOfAttackGeom()*area;
        cl += elements_[i].liftCoefficient()*area;
        cd += elements_[i].dragCoefficient()*area;
        cm += elements_[i].momentCoefficient()*area;
    }

    x /= nElements_; y /= nElements_; z /= nElements_;
    relVelMag /= totalArea;
    alphaDeg /= totalArea;
    alphaGeom /= totalArea;
    cl /= totalArea; cd /= totalArea; cm /= totalArea;

    // write time,x,y,z,rel_vel_mag,alpha_deg,alpha_geom_deg,cl,cd,cm
    *outputFile_<< time << "," << x << "," << y << "," << z << "," << relVelMag
                << "," << alphaDeg << "," << alphaGeom << "," << cl << ","
                << cd << "," << cm << endl;
}


void Foam::fv::actuatorLineSource::calcEndEffects()
{
    if (debug)
    {
        Info<< "Calculating end effects for " << name_ << endl;
    }

    scalar pi = Foam::constant::mathematical::pi;
    List<scalar> c(nElements_, 1.0); // Chord lengths
    List<scalar> alpha(nElements_, 0.1); // Geometric AoA in radians
    List<scalar> theta(nElements_); // Span distance rescaled on [0, pi]
    List<scalar> relVelMag(nElements_, 1.0);
    simpleMatrix<scalar> D(nElements_, 0.0, 0.1);
    List<scalar> A(nElements_); // Fourier coefficients
    List<scalar> circulation(nElements_);
    List<scalar> cl(nElements_);

    // Create lists from element parameters
    forAll(elements_, n)
    {
        theta[n] = elements_[n].rootDistance()*pi;
        c[n] = elements_[n].chordLength();
        //~ alpha[n] = Foam::degToRad(elements_[n].angleOfAttackGeom());
        //~ relVelMag[n] = mag(elements_[n].relativeVelocityGeom());
    }

    // Create D matrix
    forAll(elements_, i)
    {
        scalar n = i + 1;
        forAll(elements_, m)
        {
            D[m][i] = 2.0*totalLength_/(pi*c[m])*sin(n*theta[m])
                    + n*sin(n*theta[m]) / sin(theta[m]);
        }
        D.source()[i] = alpha[i];
    }
    A = D.solve();

    forAll(elements_, m)
    {
        scalar sumA = 0.0;
        forAll(elements_, i)
        {
            scalar n = i + 1;
            sumA += A[i]*sin(n*theta[m]);
        }
        circulation[m] = 2*totalLength_*relVelMag[m]*sumA;
        cl[m] = circulation[m]/(0.5*c[m]*relVelMag[m]);
    }

    // Set endEffectFactor for all elements
    List<scalar> factors = cl/Foam::max(cl);
    forAll(elements_, i)
    {
        elements_[i].setEndEffectFactor(factors[i]);
    }

    if (debug == 2)
    {
        Info<< "Debug output from actuatorLineSource::calcEndEffects:" << endl;
        Info<< "theta: " << theta << endl;
        Info<< "A: " << A << endl;
        Info<< "c: " << c << endl;
        Info<< "D.source: " << D.source() << endl;
        Info<< "D: " << D << endl;
        Info<< "cl: " << cl << endl;
        Info<< "factors:" << factors << endl;
    }
}


void Foam::fv::actuatorLineSource::harmonicPitching()
{
    // Pitch the actuator line if time has changed
    scalar t = mesh_.time().value();
    if (t != lastMotionTime_)
    {
        scalar omega = reducedFreq_*2*mag(freeStreamVelocity_)/chordLength_;
        scalar dt = mesh_.time().deltaT().value();
        scalar deltaPitch = degToRad(pitchAmplitude_)*(sin(omega*t)
                          - sin(omega*(t - dt)));
        pitch(deltaPitch);
        lastMotionTime_ = t;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::actuatorLineSource::actuatorLineSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    cellSetOption(name, modelType, dict, mesh),
    force_(vector::zero),
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
            dimForce/dimVolume,
            vector::zero
        )
    ),
    writePerf_(coeffs_.lookupOrDefault("writePerf", false)),
    lastMotionTime_(mesh.time().value()),
    endEffectsActive_(false)
{
    read(dict_);
    createElements();
    if (writePerf_)
    {
        createOutputFile();
    }
    if (forceField_.writeOpt() == IOobject::AUTO_WRITE)
    {
        forceField_.write();
    }
    // Calculate end effects
    if (endEffectsActive_)
    {
        calcEndEffects();
    }
}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

Foam::fv::actuatorLineSource::~actuatorLineSource()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::actuatorLineSource::printCoeffs() const
{
    // Print turbine properties
    Info<< "Actuator line properties:" << endl;
    Info<< "Profile data:" << endl;
    Info<< profileData_ << endl;
    Info<< "First item of element geometry:" << endl;
    Info<< elementGeometry_[0] << endl;
}


void Foam::fv::actuatorLineSource::rotate
(
    vector rotationPoint,
    vector axis,
    scalar radians
)
{
    forAll(elements_, i)
    {
        elements_[i].rotate(rotationPoint, axis, radians, true);
    }
}


void Foam::fv::actuatorLineSource::pitch(scalar radians)
{
    forAll(elements_, i)
    {
        elements_[i].pitch(radians);
    }
}


void Foam::fv::actuatorLineSource::pitch(scalar radians, scalar chordFraction)
{
    forAll(elements_, i)
    {
        elements_[i].pitch(radians, chordFraction);
    }
}


void Foam::fv::actuatorLineSource::translate(vector translationVector)
{
    forAll(elements_, i)
    {
        elements_[i].translate(translationVector);
    }
}


void Foam::fv::actuatorLineSource::setSpeed
(
    vector point,
    vector axis,
    scalar omega
)
{
    forAll(elements_, i)
    {
        elements_[i].setSpeed(point, axis, omega);
    }
}


void Foam::fv::actuatorLineSource::scaleVelocity(scalar scale)
{
    forAll(elements_, i)
    {
        elements_[i].scaleVelocity(scale);
    }
}


void Foam::fv::actuatorLineSource::setOmega(scalar omega)
{
    forAll(elements_, i)
    {
        elements_[i].setOmega(omega);
    }
}


const Foam::vector& Foam::fv::actuatorLineSource::force()
{
    return force_;
}


const Foam::volVectorField& Foam::fv::actuatorLineSource::forceField()
{
    return forceField_;
}


PtrList<Foam::fv::actuatorLineElement>& Foam::fv::actuatorLineSource::elements()
{
    return elements_;
}


Foam::vector Foam::fv::actuatorLineSource::moment(vector point)
{
    vector moment(vector::zero);
    forAll(elements_, i)
    {
        moment += elements_[i].moment(point);
    }

    if (debug)
    {
        Info<< "Moment on " << name_ << " about " << point << ": " << moment
            << endl;
    }

    return moment;
}


void Foam::fv::actuatorLineSource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    // If harmonic pitching is active, do harmonic pitching
    if (harmonicPitchingActive_)
    {
        harmonicPitching();
    }

    // Zero out force field
    forceField_ *= dimensionedScalar("zero", forceField_.dimensions(), 0.0);

    // Zero the total force vector
    force_ = vector::zero;

    forAll(elements_, i)
    {
        elements_[i].addSup(eqn, forceField_);
        force_ += elements_[i].force();
    }

    Info<< "Force (per unit density) on " << name_ << ": "
        << endl << force_ << endl << endl;

    // Check dimensions on force field and correct if necessary
    if (forceField_.dimensions() != eqn.dimensions()/dimVolume)
    {
        forceField_.dimensions().reset(eqn.dimensions()/dimVolume);
    }

    // Add source to eqn
    eqn += forceField_;

    // Write performance to file
    if (writePerf_ and Pstream::master())
    {
        writePerf();
    }
}


void Foam::fv::actuatorLineSource::addSup
(
    fvMatrix<scalar>& eqn,
    const label fieldI
)
{
    // If harmonic pitching is active, do harmonic pitching
    if (harmonicPitchingActive_)
    {
        harmonicPitching();
    }

    const volVectorField& U = mesh_.lookupObject<volVectorField>("U");

    word fieldName = fieldNames_[fieldI];

    Info<< endl << "Adding " << fieldName << " from " << name_ << endl << endl;
    forAll(elements_, i)
    {
        elements_[i].calculateForce(U);
        elements_[i].addTurbulence(eqn, fieldName);
    }
}


void Foam::fv::actuatorLineSource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    // If harmonic pitching is active, do harmonic pitching
    if (harmonicPitchingActive_)
    {
        harmonicPitching();
    }

    // Zero out force field
    forceField_ *= dimensionedScalar("zero", forceField_.dimensions(), 0.0);

    // Zero the total force vector
    force_ = vector::zero;

    forAll(elements_, i)
    {
        elements_[i].addSup(rho, eqn, forceField_);
        force_ += elements_[i].force();
    }

    Info<< "Force on " << name_ << ": " << endl << force_ << endl << endl;

    // Check dimensions of force field and correct if necessary
    if (forceField_.dimensions() != eqn.dimensions()/dimVolume)
    {
        forceField_.dimensions().reset(eqn.dimensions()/dimVolume);
    }

    // Add source to eqn
    eqn += forceField_;

    // Write performance to file
    if (writePerf_ and Pstream::master())
    {
        writePerf();
    }
}


// ************************************************************************* //
