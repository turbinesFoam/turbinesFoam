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
        coeffs_.lookup("hubHeight") >> hubHeight_;
        freeStreamDirection_ = freeStreamVelocity_/mag(freeStreamVelocity_);
        endEffectsActive_ = coeffs_.lookupOrDefault("endEffects", false);

        // Multi-phase parameters, if present
        multiPhase_ = coeffs_.lookupOrDefault("multiPhase", false);
        if(multiPhase_)
        {
            coeffs_.lookup("phaseName") >> phaseName_;
        }

        // Read harmonic pitching parameters if present
        dictionary pitchDict = coeffs_.subOrEmptyDict("harmonicPitching");
        harmonicPitchingActive_ = pitchDict.lookupOrDefault("active", false);
        reducedFreq_ = pitchDict.lookupOrDefault("reducedFreq", 0.0);
        pitchAmplitude_ = pitchDict.lookupOrDefault("amplitude", 0.0);

        //New- Read wind pitching parameters if present
        dictionary pitchDict = coeffs_.subOrEmptyDict("windPitching");
        windPitchingActive_ = pitchDict.lookupOrDefault("active", false);

        // Read prescribed motion dictionary if present
        prescribedMotionDict_ = coeffs_.subOrEmptyDict("prescribedMotion");
        prescribedMotionActive_ = prescribedMotionDict_.lookupOrDefault("active", false);
        if (prescribedMotionActive_)
        {
            prescribedMotionInitialise();
        }

        // Read rigid body motion parameters if present
        rigidBodyMotionDict_ = coeffs_.subOrEmptyDict("rigidBodyMotion");
        rigidBodyMotionActive_ = rigidBodyMotionDict_.lookupOrDefault("active", false);
        rigidBodyAligned_ = rigidBodyMotionDict_.lookupOrDefault("isAligned", false);
        if(rigidBodyMotionActive_)
        {
            refOrientation_ =  rigidBodyMotionDict_.lookupOrDefault("principalOrientation", tensor::I);
            if(equal(det(refOrientation_), 0.0))
            {
                Info << "Principal orientation must have non-zero determinant." << endl;
                Info << "Rigid body motion is disabled." << endl;
                rigidBodyMotionActive_ = false;
            }
        }

        // Check whether rigid body and prescribed harmonic are both active
        if(rigidBodyMotionActive_ && prescribedMotionActive_)
        {
            Info<< "Both 'prescribedMotion' and 'rigidBodyMotion' are set to active. " <<
            "Since only one type of motion is allowed, they will be disabled." << endl;
            rigidBodyMotionActive_ =   false;
            prescribedMotionActive_  =   false;
        }


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

    *outputFile_<< "time,x,y,z,rel_vel_mag,motion_vel_mag,alpha_deg,alpha_geom_deg,cl,cd,cm"
                << endl;
}


void Foam::fv::actuatorLineSource::createElements()
{
    elements_.setSize(nElements_);

    label nGeometryPoints = elementGeometry_.size();
    label nGeometrySegments = nGeometryPoints - 1;
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
        if (i > 0)
        {
            spanLengths[i - 1] = mag(points[i] - points[i-1]);
            totalLength_ += spanLengths[i - 1];
        }
        // Read span direction
        x = elementGeometry_[i][1][0];
        y = elementGeometry_[i][1][1];
        z = elementGeometry_[i][1][2];
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

    // Store blade root and tip locations for distance calculations
    vector rootLocation = points[0];
    vector tipLocation = points[nGeometryPoints - 1];

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
        vector point1 = points[geometrySegmentIndex];
        vector point2 = points[geometrySegmentIndex + 1];
        vector segment = point2 - point1;
        position = point1
                 + segment/nElementsPerSegment*pointIndex
                 + segment/nElementsPerSegment/2;

        // Linearly interpolate chordLength
        scalar chordLength1 = chordLengths[geometrySegmentIndex];
        scalar chordLength2 = chordLengths[geometrySegmentIndex + 1];
        scalar deltaChordTotal = chordLength2 - chordLength1;
        chordLength = chordLength1
                    + deltaChordTotal/nElementsPerSegment*pointIndex
                    + deltaChordTotal/nElementsPerSegment/2;

        // Linearly interpolate spanDirection
        vector spanDir1 = spanDirs[geometrySegmentIndex];
        vector spanDir2 = spanDirs[geometrySegmentIndex + 1];
        vector deltaSpanTotal = spanDir2 - spanDir1;
        spanDirection = spanDir1
                      + deltaSpanTotal/nElementsPerSegment*pointIndex
                      + deltaSpanTotal/nElementsPerSegment/2;

        // Linearly interpolate section pitch
        scalar pitch1 = pitches[geometrySegmentIndex];
        scalar pitch2 = pitches[geometrySegmentIndex + 1];
        scalar deltaPitchTotal = pitch2 - pitch1;
        pitch = pitch1
              + deltaPitchTotal/nElementsPerSegment*pointIndex
              + deltaPitchTotal/nElementsPerSegment/2;

        // Linearly interpolate chord mount
        scalar cm1 = chordMounts[geometrySegmentIndex];
        scalar cm2 = chordMounts[geometrySegmentIndex + 1];
        scalar deltaCmTotal = cm2 - cm1;
        chordMount = cm1 + deltaCmTotal/nElementsPerSegment*pointIndex
                   + deltaCmTotal/nElementsPerSegment/2;

        // Linearly interpolate element velocity
        vector vel1 = initialVelocities[geometrySegmentIndex];
        vector vel2 = initialVelocities[geometrySegmentIndex + 1];
        vector deltaVelTotal = vel2 - vel1;
        initialVelocity = vel1
                        + deltaVelTotal/nElementsPerSegment*pointIndex
                        + deltaVelTotal/nElementsPerSegment/2;

        // Linearly interpolate chordDirection
        vector chordDir1 = chordRefDirs[geometrySegmentIndex];
        vector chordDir2 = chordRefDirs[geometrySegmentIndex + 1];
        vector deltaChordDirTotal = chordDir2 - chordDir1;
        chordDirection = chordDir1
                       + deltaChordDirTotal/nElementsPerSegment*pointIndex
                       + deltaChordDirTotal/nElementsPerSegment/2;

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
        dict.add("multiPhase", multiPhase_);
        dict.add("phaseName", phaseName_);

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
    scalar motionVelMag = 0.0;
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
        motionVelMag += mag(elements_[i].motionVelocity())*area;
        alphaDeg += elements_[i].angleOfAttack()*area;
        alphaGeom += elements_[i].angleOfAttackGeom()*area;
        cl += elements_[i].liftCoefficient()*area;
        cd += elements_[i].dragCoefficient()*area;
        cm += elements_[i].momentCoefficient()*area;
    }

    x /= nElements_; y /= nElements_; z /= nElements_;
    relVelMag /= totalArea;
    motionVelMag /= totalArea;
    alphaDeg /= totalArea;
    alphaGeom /= totalArea;
    cl /= totalArea; cd /= totalArea; cm /= totalArea;

    // write time,x,y,z,rel_vel_mag,motion_vel_mag,alpha_deg,alpha_geom_deg,cl,cd,cm
    *outputFile_<< time << "," << x << "," << y << "," << z << "," << relVelMag
                << "," << motionVelMag << ","  << alphaDeg << "," << alphaGeom 
                << "," << cl << "," << cd << "," << cm << endl;
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
        pitch(deltaPitch, 0);
        lastMotionTime_ = t;
    }
}

void Foam::fv::actuatorLineSource::windPitching(const fvMesh& mesh)
{
    // Pitch the actuator line if wind speed reaches threshold and time has changed
    scalar t = mesh_.time().value();
    if (t != lastMotionTime_)
    {   
        readRigidBodyDict(mesh);
        rigidBodyDict_.lookup("velocity") >> linearVelocity_;
        rigidBodyDict_.lookup("omega") >> angularVelocity_;
    
        // Calculate the relative speed to determine whether to pitch or not
        scalar Uhub = freeStreamVelocity_.x();
        scalar Umotionx = linearVelocity_.x();  //The speed brought by the Surge
        scalar Upitch = angularVelocity_.y();  //The relative velocity caused by pitch motion at the hub
        Urel_ = Uhub - Umotionx - Upitch*hubHeight_; 
        
        scalar targetPitch = 0;

        // Calculate the target pitch angle based on the relative velocity
        if (Urel_ >= 10.673 &&  Urel_ < 11.17)
        {
            targetPitch = 0.512;
        }
        if (Urel_ >= 11.17 &&  Urel_ < 11.7)
        {
            targetPitch = 3.72;
        }
        if (Urel_ >= 11.7 &&  Urel_ < 12.26)
        {
            targetPitch = 5.39;
        }
        if (Urel_ >= 12.26 &&  Urel_ < 12.84)
        {
            targetPitch = 6.76;
        }
        if (Urel_ >= 12.84 &&  Urel_ < 13.46)
        {
            targetPitch = 7.98;
        }
        if (Urel_ >= 13.46 &&  Urel_ < 14.10)
        {
            targetPitch = 9.12;
        }
        if (Urel_ >= 14.10 &&  Urel_ < 14.77)
        {
            targetPitch = 10.2;
        }
        if (Urel_ >= 14.77 &&  Urel_ < 15.47)
        {
            targetPitch = 11.23;
        }
        if (Urel_ >= 15.47 &&  Urel_ < 16.18)
        {
            targetPitch = 12.23;
        }
        if (Urel_ >= 16.18 &&  Urel_ < 16.92)
        {
            targetPitch = 13.21;
        }
        if (Urel_ >= 16.92 &&  Urel_ < 17.67)
        {
            targetPitch = 14.16;
        }
        if (Urel_ >= 17.67 &&  Urel_ < 18.44)
        {
            targetPitch = 15.10;
        }
        if (Urel_ >= 18.44 &&  Urel_ < 19.23)
        {
            targetPitch = 16.02;
        }
        if (Urel_ >= 19.23 &&  Urel_ < 20.03)
        {
            targetPitch = 16.93;
        }
        if (Urel_ >= 20.03 &&  Urel_ < 20.84)
        {
            targetPitch = 17.82;
        }
        if (Urel_ >= 20.84 &&  Urel_ < 21.66)
        {
            targetPitch = 18.70;
        }
        if (Urel_ >= 21.66 &&  Urel_ < 22.48)
        {
            targetPitch = 19.56;
        }
        if (Urel_ >= 22.48 &&  Urel_ < 23.32)
        {
            targetPitch = 20.41;
        }
        if (Urel_ >= 23.32 &&  Urel_ < 24.16)
        {
            targetPitch = 21.25;
        }
        if (Urel_ >= 24.16 &&  Urel_ < 25)
        {
            targetPitch = 22.07;
        }
        if (Urel_ >= 25)
        {
            targetPitch = 22.88;
        }

        scalar deltaPitch = targetPitch - lastPitchAngle_;

        // Apply pitch
        pitch(deltaPitch);

        //log pitch
        lastPitchAngle_ = targetPitch;
        lastMotionTime_ = t;
    }
}

void Foam::fv::actuatorLineSource::prescribedMotionInitialise()
{
    prescribedMotionDict_.lookup("motionType") >> prescribedMotionType_;
    prescribedMotionCoeffs_ = prescribedMotionDict_.subOrEmptyDict
    (
        prescribedMotionType_ + "Coeffs"
    );

    if (prescribedMotionType_ == "harmonic1st" || prescribedMotionType_ == "customMotion")
    {
        //- Center of rotation of the prescribed/rigidBody motion
        rotCenter_ = prescribedMotionDict_.lookupOrDefault("rotationCenter", vector::zero);
        
        //- Orientation at t=0 of the prescribed/rigidBody motion
        rot0_ = prescribedMotionDict_.lookupOrDefault("initialRotation", vector::zero); 
        rot0_ *= M_PI / 180; // Transform deg into rad  

        //Initialise orientation according to prescribed motion
        orientation_ =      rotAngles2Matrix(rot0_);
        prevOrientation_ =  tensor::I;
        prevRotCenter_ =    rotCenter_;

        if (debug)
        {
            Info << "Chosen '"<< prescribedMotionType_ << "' prescribed motion type for " << name_ << "." << endl;
        }
    }
    else
    {
        Info << "Specified prescribed motion type '" << prescribedMotionType_ << "' not implemented. "<< endl;
        Info << "Disabling prescribed motion. "<< endl;
        prescribedMotionActive_ = false;
    }
}

void Foam::fv::actuatorLineSource::prescribeMotion()
{
    scalar t = mesh_.time().value();
    // If time has changed, move the actuator line according to prescribed motion
    if (t != lastMotionTime_)
    {
        Info << "Applying prescribed motion at time t =  "<< t << endl;
        scalar dt = mesh_.time().deltaT().value();

        // Update previous rotation.
        prevOrientation_ = orientation_;

        // Compute prescribed motion parameters
        // translation: change in position from the previous time-step: r(t)-r(t-dt)
        // rotation: euler angles defining the absolute orientation 
        // velocity: linear velocity
        // omega: angular velocity

        vector translation(vector::zero), rotation(vector::zero), velocity(vector::zero), omega(vector::zero);
        prescribedMotionParameters(translation, rotation, velocity, omega, t, dt);

        // Rotation matrix: from floater (F) to inertial (I) frame
        orientation_ = rotAngles2Matrix(rotation);
        
        // omega is given in floater (F) frame and thus
        // has to be transformed into inertial (I)
        omega = orientation_ & omega;

        // Move the actuator line according to the computed values
        applyMotion(translation, velocity, omega);
        lastMotionTime_ = t;
    }
}

void Foam::fv::actuatorLineSource::prescribedMotionParameters
(
    vector &translation,
    vector &rotation,
    vector &velocity,
    vector &omega,
    const scalar t,
    const scalar dt
)
{
    // First-order harmonic motion
    if (prescribedMotionType_ == "harmonic1st")
    {
        //- Amplitude of prescribed translation motion (x y z, inertial frame) in meters 
        vector harmonicTraAmp =   prescribedMotionCoeffs_.lookupOrDefault("translationAmplitude", vector::zero);
        //- Frequency of prescribed translation motion in Hz
        vector harmonicTraFreq =  prescribedMotionCoeffs_.lookupOrDefault("translationFrequency", vector::zero);
        //- Amplitude of prescribed rotation motion (roll pitch yaw, rotating frame) in deg
        vector harmonicRotAmp =   prescribedMotionCoeffs_.lookupOrDefault("rotationAmplitude", vector::zero);
        //- Frequency of prescribed rotation motion in Hz
        vector harmonicRotFreq =  prescribedMotionCoeffs_.lookupOrDefault("rotationFrequency", vector::zero);
        
        harmonicTraFreq *= 2*M_PI; // Transform Hz into rad/s
        harmonicRotAmp *= M_PI / 180; // Transform deg into rad     
        harmonicRotFreq *= 2*M_PI; // Transform Hz into rad/s      

        forAll(translation, i)
        {
            translation[i] = harmonicTraAmp[i] * (sin(harmonicTraFreq[i]*t)
                       - sin(harmonicTraFreq[i]*(t-dt))) ;                         
            rotation[i] = rot0_[i] + harmonicRotAmp[i] * sin(harmonicRotFreq[i]*t);                  
            velocity[i] = harmonicTraFreq[i] * harmonicTraAmp[i]
                           * cos(harmonicTraFreq[i] * t);
            omega[i] = harmonicRotFreq[i] * harmonicRotAmp[i]
                           * cos(harmonicRotFreq[i] * t);
        }
    }
    else if (prescribedMotionType_ == "customMotion")
    {
        // Custom motion, define the desired function below
        // The case: prescribed motion in pitch and surge only
        scalar w = M_PI/4; // Frequency

        scalar surge_p = M_PI/2 + 10.4776*M_PI/180; // Surge phase
        scalar surge_a = 0.0430; // Surge velocity amplitude
        scalar surge_0 = -0.0100; // Surge-zero
        scalar vel_0 = -0.006699; // Surge velocity - zero

        scalar surge_vel = vel_0 + surge_a*sin(w*t + surge_p);
        scalar surge = surge_0 + vel_0*t - (surge_a/w)*cos(w*t+surge_p);
        scalar surge_prev = surge_0 + vel_0*(t-dt) - (surge_a/w)*cos(w*(t-dt)+surge_p);

        scalar pitch_p = M_PI/2 + 140.8746*M_PI/180; // Pitch phase
        scalar pitch_a = 0.1722; // Pitch velocity amplitude
        scalar pitch_0 = 2.218; // Pitch-zero

        scalar pitch = pitch_0 + pitch_a*sin(w*t + pitch_p);
        scalar pitch_vel = w*pitch_a*cos(w*t + pitch_p);

        pitch *= M_PI / 180; // Transform deg into rad     
        pitch_vel *= M_PI / 180; // Transform deg/s into rad/s 

        translation[0] = surge - surge_prev;
        velocity[0] = surge_vel;
        rotation[1] = pitch;
        omega[1] = pitch_vel;
    }
}

void Foam::fv::actuatorLineSource::readRigidBodyDict(const fvMesh& mesh)
{
    // If rigid body dictionary exists
    if(mesh.time().foundObject<IOdictionary>("sixDoFMotion"))
    {
        // Access dictionary
        const dictionary& motionDict = mesh.time().lookupObject<IOdictionary>("sixDoFMotion");
        // Save its reference
        rigidBodyDict_ = motionDict;
        if (debug)
        {
            Info << "Chosen 'rigidBody' motion type for " << name_ << "." << endl;
        };
    }
    else
    {
        rigidBodyMotionActive_ = false;
        Info << "Rigid body IOdictionary ('sixDoFMotion') could not be accessed." << endl;
        Info << "Motion based on rigid body is thus disabled." << endl;
    }
}


void Foam::fv::actuatorLineSource::rigidBodyInitialise(const fvMesh& mesh)
{    
    if(rigidBodyMotionActive_)
    {
        readRigidBodyDict(mesh);
        // If 'sixDoFMotion' IOdictionary has been accessed
        if(rigidBodyMotionActive_)
        {
            // Access its content
            rigidBodyDict_.lookup("centreOfRotation") >> rotCenter_;
            rigidBodyDict_.lookup("orientation") >> rBOrientation_;
            orientation_ = rBOrientation_ & refOrientation_;
            prevOrientation_ = tensor::I;
            prevRotCenter_ = rotCenter_;
        }  
    }    
}


void Foam::fv::actuatorLineSource::rigidBodyMotion(const fvMesh& mesh)
{
    // Update previous values
    prevOrientation_ = orientation_;
    prevRotCenter_ = rotCenter_;

    // Read the IOdictionary
    readRigidBodyDict(mesh);

    // Update current rotation center
    rigidBodyDict_.lookup("centreOfRotation") >> rotCenter_;

    // Rigid body orientation 'Q' maps from local to global.
    // This is the same as going from floater (F) to inertial
    // (I) frames. Thus, no need to transpose it. 
    rigidBodyDict_.lookup("orientation") >> rBOrientation_;

    // Rigid body orientation 'Q' is given wrt the orientation
    // aligned with inertia principal axes (refOrientation).
    // So, first apply refOrientation and then Q.
    orientation_ = rBOrientation_ & refOrientation_;

    // Get velocities
    vector translation, velocity, omega;
    translation = rotCenter_ - prevRotCenter_;
    rigidBodyDict_.lookup("velocity") >> velocity;
    rigidBodyDict_.lookup("omega") >> omega;

    // Move the actuator line according to the computed values
    applyMotion(translation, velocity, omega);

    if(debug)
    {
        Info << "Rigid body motion." << endl;
        Info << "Rigid body orientation: " << rBOrientation_ << endl;
        Info << "Overall orientation: " << orientation_ << endl;
        Info << "Rotation centre: " << rotCenter_ << endl;
        Info << "Velocity: " << velocity << endl;
        Info << "Omega: " << omega << endl;
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
    pitchRelativeVelocity_(vector::zero),
    freeStreamVelocity_(vector::zero),
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
    rigidBodyInitialise(mesh);
    createElements();
    motionAlign();
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
    lastPitchAngle_ = 0.0;        // 初始桨距角为 0
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


void Foam::fv::actuatorLineSource::rotate
(
    const vector &rotationPoint,
    const tensor &rotMatrix
)
{
    forAll(elements_, i)
    {
        elements_[i].rotate(rotationPoint, rotMatrix, true);
    }
}


void Foam::fv::actuatorLineSource::pitch(scalar radians)
{
    forAll(elements_, i)
    {
        elements_[i].pitch(radians);
    }
}

const Foam::vector& Foam::fv::actuatorLineSource::pitchRelativeVelocity()
{
    return pitchRelativeVelocity_;
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


void Foam::fv::actuatorLineSource::motionRotate()
{
    // Total rotation matrix: 
    // go back to unrotated orientation 
    // and rotate again with rotMatrix
    // R = Rn * R(n-1)^T
    tensor totalRotMatrix = orientation_ & prevOrientation_.T();

    // Execute rotation routine only if rotation
    // angle is not zero. From the formula:
    // tr(R) = 1 + 2cos(angle)
    // angle==0 if tr(R) ==3
    if(!equal(tr(totalRotMatrix), scalar(3.0)))
    {
        // Perform rotation wrt the rotation center
        // from the previous timestep
        // R * ( x_(n-1) - c_(n-1) )
        rotate(rotCenter_, totalRotMatrix);
        if(debug)
        {
            Info<< "Accounting for prescribed/rigidBody motion... " << endl;
            Info<< "Actuator line rotation matrix: " << orientation_ << endl;
            Info<< "Rotation center: " << rotCenter_ << endl;
        }
    }  
}


void Foam::fv::actuatorLineSource::applyMotion
(
    const vector &translation,
    const vector &velocity,
    const vector &omega
)
{
    // Motion parameters
    // translation: change in position from the previous time-step: r(t)-r(t-dt)
    // velocity: linear velocity
    // omega: angular velocity

    // Rotate according to prescribed/rigidBody motion
    motionRotate();

    // Now translate 
    translate(translation);
    
    // Update rotation center
    prevRotCenter_ = rotCenter_;
    rotCenter_ += translation;

    // Set velocity according to prescribed/rigidBody motion
    setMotionVelocity(velocity);

    if(debug)
    {
        Info<< "Translating: " << translation << " [m]" << endl; 
    }

    // Add velocity due to omega
    scalar omegaNorm = mag(omega);
    if (omegaNorm > SMALL)
    {
        // Rotation vector is given by the angular speed
        vector omegaAxis = omega / omegaNorm;
        addMotionOmega(rotCenter_, omegaAxis, omegaNorm);
    }
}


void Foam::fv::actuatorLineSource::motionAlign()
{
    if(prescribedMotionActive_ || (rigidBodyMotionActive_ && !rigidBodyAligned_))
    {
        // Rotate according to prescribed/rigidBody motion at t=0
        motionRotate();
    }
}


tensor Foam::fv::actuatorLineSource::rotAngles2Matrix(const vector &rotAngles)
{
    scalar c1 = cos(rotAngles.z());
    scalar s1 = sin(rotAngles.z());
    scalar c2 = cos(rotAngles.y());
    scalar s2 = sin(rotAngles.y());
    scalar c3 = cos(rotAngles.x());
    scalar s3 = sin(rotAngles.x());

    // Compute rotation matrix as the intrinsic rotation (z1 y2 x3) 
    // with Tait-Bryan angles roll, pitch and yaw
    tensor rotMatrix;
    rotMatrix.xx() = c1*c2;
    rotMatrix.xy() = c1*s2*s3 - c3*s1;
    rotMatrix.xz() = s1*s3 + c1*c3*s2;
    rotMatrix.yx() = c2*s1;
    rotMatrix.yy() = c1*c3 + s1*s2*s3;
    rotMatrix.yz() = c3*s1*s2 - c1*s3;
    rotMatrix.zx() = -s2;
    rotMatrix.zy() = c2*s3;
    rotMatrix.zz() = c2*c3;

    return rotMatrix;
}

void Foam::fv::actuatorLineSource::rotMatrix2Axis
(
    const tensor &rotMatrix,
    vector &axis,
    scalar &angle
)
{
    // Transform the rotation matrix into an
    // axis-angle rotation using Rodrigues formula.
    // Algorithm taken from 
    // https://courses.cs.duke.edu/fall13/compsci527/notes/rodrigues.pdf

    tensor A = skew(rotMatrix);
    vector rho = vector(A.zy(), A.xz(), A.yx());
    scalar s = mag(rho);
    scalar c = (tr(rotMatrix) - 1) / 2;

    // In case sin(angle) = 0 and cos(angle) = 1,
    // there is no rotation
    if (s==0 && c==1)
    {
        // Dummy axis, rotation angle is zero
        axis = vector(1, 0, 0);
        angle = 0.0;
    }
    // In case sin(angle) = 0 and cos(angle) = -1,
    // angle could be +-pi. 
    else if (s==0 && c==-1)
    {
        // Pick a non-zero column of V = R-I.
        // We know that at least one column
        // must have non-zero norm
        tensor V = rotMatrix - tensor::I;
        vector v = V.cx();
        if (mag(v) < SMALL)
        {
            v = V.cy();
        }
        if (mag(v) < SMALL)
        {
            v = V.cz();
        } 
        axis = v/mag(v);

        // To ensure angle uniqueness, distinguish
        // between +-pi cases. This might not be
        // necessary for the present application. 
        scalar r1 = axis.x();
        scalar r2 = axis.y();
        scalar r3 = axis.z();
        if ((mag(r1)<SMALL && mag(r2)<SMALL && r3<0)
            || (mag(r1)<SMALL && r2<0) || (r1<0))
        {
            angle = -M_PI;
        }
        else
        {
            angle = M_PI;
        }
    }
    // General case, where sin(angle) =/ 0
    else
    {
        axis = rho / s;
        angle = atan2(s,c);
    }

}


tensor Foam::fv::actuatorLineSource::rotAxis2Matrix
(
    const vector &axis,
    const scalar &angle
)
{
    // Declare and define the rotation matrix (from SOWFA)
    tensor RM;
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

    return RM;
}


void Foam::fv::actuatorLineSource::setMotionVelocity(const vector &velocityVector)
{
    forAll(elements_, i)
    {
        elements_[i].setMotionVelocity(velocityVector);
    }
}


void Foam::fv::actuatorLineSource::addMotionVelocity(const vector &velocityVector)
{
    forAll(elements_, i)
    {
        elements_[i].addMotionVelocity(velocityVector);
    }
}


void Foam::fv::actuatorLineSource::addMotionOmega
(
    const vector &point,
    const vector &axis,
    const scalar &omega
)
{
    forAll(elements_, i)
    {
        elements_[i].addMotionOmega(point, axis, omega);
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


const Foam::tensor& Foam::fv::actuatorLineSource::orientation()
{
    return orientation_;
}


const Foam::tensor& Foam::fv::actuatorLineSource::prevOrientation()
{
    return prevOrientation_;
}


const Foam::vector& Foam::fv::actuatorLineSource::rotCenter()
{
    return rotCenter_;
}


const Foam::vector& Foam::fv::actuatorLineSource::prevRotCenter()
{
    return prevRotCenter_;
}


const Foam::vector Foam::fv::actuatorLineSource::motionTranslation()
{
    vector translation = rotCenter_ - prevRotCenter_;
    return translation;
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


void Foam::fv::actuatorLineSource::addSup //- Source term to momentum equation
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

    // If wind pitching is active, do wind pitching
    if (windPitchingActive_)
    {
        windPitching();
    }

    // If prescribed motion is active, move the actuator line accordingly
    if (prescribedMotionActive_)
    {
        prescribeMotion();
    }

    // If rigid body motion is active, move the actuator line accordingly
    if (rigidBodyMotionActive_)
    {
        rigidBodyMotion(mesh_);
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


void Foam::fv::actuatorLineSource::addSup //- Source term to turbulence scalars
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

    if (windPitchingActive_)
    {
        windPitching();
    }

    // If prescribed motion is active, move the actuator line accordingly
    if (prescribedMotionActive_)
    {
        prescribeMotion();
    }

    // If rigid body motion is active, move the actuator line accordingly
    if (rigidBodyMotionActive_)
    {
        rigidBodyMotion(mesh_);
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


void Foam::fv::actuatorLineSource::addSup //- Source term to compressible momentum equation
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

    if (windPitchingActive_)
    {
        windPitching();
    }

    // If prescribed motion is active, move the actuator line accordingly
    if (prescribedMotionActive_)
    {
        prescribeMotion();
    }

    // If rigid body motion is active, move the actuator line accordingly
    if (rigidBodyMotionActive_)
    {
        rigidBodyMotion(mesh_);
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
