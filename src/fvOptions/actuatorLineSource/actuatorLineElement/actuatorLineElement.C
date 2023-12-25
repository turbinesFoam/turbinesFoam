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

#include "actuatorLineElement.H"
#include "addToRunTimeSelectionTable.H"
#include "geometricOneField.H"
#include "fvMatrices.H"
#include "syncTools.H"
#include "unitConversion.H"

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
    dict_.lookup("chordRefDirection") >> chordRefDirection_;
    dict_.lookup("chordMount") >> chordMount_;
    dict_.lookup("spanLength") >> spanLength_;
    dict_.lookup("spanDirection") >> spanDirection_;
    dict_.lookup("freeStreamVelocity") >> freeStreamVelocity_;
    freeStreamDirection_ = freeStreamVelocity_/mag(freeStreamVelocity_);
    dict_.lookup("rootDistance") >> rootDistance_;
    dict_.lookup("velocitySampleRadius") >> velocitySampleRadius_;
    dict_.lookup("nVelocitySamples") >> nVelocitySamples_;

    // Create dynamic stall model if found
    if (dict_.found("dynamicStall"))
    {
        dictionary dsDict = dict_.subDict("dynamicStall");
        word dsName;
        dsDict.lookup("dynamicStallModel") >> dsName;
        dynamicStall_ = dynamicStallModel::New
        (
            dsDict,
            dsName,
            mesh_.time(),
            profileData_
        );
        dsDict.lookup("active") >> dynamicStallActive_;
    }

    // Read flow curvature correction subdictionary
    if (dict_.found("flowCurvature"))
    {
        dictionary fcDict = dict_.subDict("flowCurvature");
        flowCurvatureActive_ = fcDict.lookupOrDefault("active", false);
        word defaultName = "none";
        flowCurvatureModelName_ = fcDict.lookupOrDefault
        (
            "flowCurvatureModel",
            defaultName
        );
    }

    // Read nu from object registry
    const dictionary& transportProperties = mesh_.lookupObject<IOdictionary>
    (
        "transportProperties"
    );
    dimensionedScalar nu;
    transportProperties.lookup("nu") >> nu;
    nu_ = nu.value();

    // Read writePerf switch
    dict_.lookup("writePerf") >> writePerf_;

    if (debug)
    {
        Info<< "actuatorLineElement properties:" << endl;
        Info<< "Position: " << position_ << endl;
        Info<< "rootDistance: " << rootDistance_ << endl;
        Info<< "chordLength: " << chordLength_ << endl;
        Info<< "chordDirection: " << chordDirection_ << endl;
        Info<< "spanLength: " << spanLength_ << endl;
        Info<< "spanDirection: " << spanDirection_ << endl;
        Info<< "writePerf: " << writePerf_ << endl;
    }
}


void Foam::fv::actuatorLineElement::rotateVector
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


Foam::label Foam::fv::actuatorLineElement::findCell
(
    const point& location
)
{
    if (Pstream::parRun())
    {
        if (meshBoundBox_.containsInside(location))
        {
            if (debug)
            {
                Pout<< "Looking for cell containing " << location
                    << " inside bounding box:" << endl
                    << meshBoundBox_ << endl;
            }
            return mesh_.findCell(location);
        }
        else
        {
            if (debug)
            {
                Pout<< "Cell not inside " << meshBoundBox_ << endl;
            }
            return -1;
        }
    }
    else
    {
        return mesh_.findCell(location);;
    }
}


void Foam::fv::actuatorLineElement::lookupCoefficients()
{
    liftCoefficient_ = profileData_.liftCoefficient(angleOfAttack_);
    dragCoefficient_ = profileData_.dragCoefficient(angleOfAttack_);
    momentCoefficient_ = profileData_.momentCoefficient(angleOfAttack_);
}


Foam::scalar Foam::fv::actuatorLineElement::calcProjectionEpsilon()
{
    // Lookup Gaussian coeffs from profileData dict if present
    dictionary GaussianCoeffs = profileData_.dict().subOrEmptyDict
    (
        "GaussianCoeffs"
    );
    scalar chordFactor = GaussianCoeffs.lookupOrDefault("chordFactor", 0.25);
    scalar dragFactor = GaussianCoeffs.lookupOrDefault("dragFactor", 1.0);
    scalar meshFactor = GaussianCoeffs.lookupOrDefault("meshFactor", 2.0);

    // Provide ideal epsilon target for lift based on chord length
    scalar epsilonLift = chordFactor*chordLength_;

    // Epsilon based on drag/momentum thickness
    scalar epsilonDrag = dragFactor*dragCoefficient_*chordLength_/2.0;

    // Threshold is based on lift or drag, whichever is larger
    scalar epsilonThreshold = Foam::max(epsilonLift, epsilonDrag);

    scalar epsilon = VGREAT;
    scalar epsilonMesh = VGREAT;
    const scalarField& V = mesh_.V();
    label posCellI = findCell(position_);
    if (posCellI >= 0)
    {
        // Projection width based on local cell size (from Troldborg (2008))
        epsilonMesh = 2.0*Foam::cbrt(V[posCellI]);
        epsilonMesh *= meshFactor; // Cell could have non-unity aspect ratio

        if (epsilonMesh > epsilonThreshold)
        {
            epsilon = epsilonMesh;
        }
        else
        {
            epsilon = epsilonThreshold;
        }
    }

    // Reduce epsilon over all processors
    reduce(epsilon, minOp<scalar>());

    // If epsilon is not reduced, position is not in the mesh
    if (not (epsilon < VGREAT))
    {
        // Raise fatal error since mesh size cannot be detected
        FatalErrorIn("void actuatorLineElement::applyForceField()")
            << "Position of " << name_  << " not found in mesh"
            << abort(FatalError);
    }

    if (debug)
    {
        reduce(epsilonMesh, minOp<scalar>());
        word epsilonMethod;
        if (epsilon == epsilonLift)
        {
            epsilonMethod = "lift-based";
        }
        else if (epsilon == epsilonDrag)
        {
            epsilonMethod = "drag-based";
        }
        else if (epsilon == epsilonMesh)
        {
            epsilonMethod = "mesh-based";
        }
        Info<< "    epsilon (" << epsilonMethod << "): " << epsilon << endl;
    }

    return epsilon;
}


void Foam::fv::actuatorLineElement::correctFlowCurvature
(
    scalar& angleOfAttackRad
)
{
    if (debug)
    {
        Info<< "    Correcting for flow curvature with "
            << flowCurvatureModelName_ << " model" << endl;
    }

    if (flowCurvatureModelName_ == "Goude")
    {
        angleOfAttackRad += omega_*chordLength_/(2*mag(relativeVelocity_));
    }
    else if (flowCurvatureModelName_ == "MandalBurton")
    {
        // Calculate relative velocity at leading and trailing edge
        vector relativeVelocityLE = inflowVelocity_ - velocityLE_;
        vector relativeVelocityTE = inflowVelocity_ - velocityTE_;

        // Calculate angle of attack at leading and trailing edge
        scalar alphaLE = asin((planformNormal_ & relativeVelocityLE)
                       / (mag(planformNormal_)*mag(relativeVelocityLE)));
        scalar alphaTE = asin((planformNormal_ & relativeVelocityTE)
                       / (mag(planformNormal_)*mag(relativeVelocityTE)));

        scalar beta = alphaTE - alphaLE;

        angleOfAttackRad += atan2((1.0 - cos(beta/2.0)), sin(beta/2.0));
    }
    else if (flowCurvatureModelName_ == "constantOffset")
    {
        dictionary fcDict = dict_.subDict("flowCurvature");
        dictionary coeffs = fcDict.subDict(flowCurvatureModelName_ + "Coeffs");
        scalar offsetDeg = 0.0;
        coeffs.lookup("offsetDeg") >> offsetDeg;
        angleOfAttackRad += degToRad(offsetDeg);
    }
}


void Foam::fv::actuatorLineElement::multiplyForceRho
(
    const volScalarField& rho
)
{
    // Lookup local density
    label cellI = findCell(position_);
    scalar localRho = VGREAT;
    if (cellI >= 0)
    {
        localRho = rho[cellI];
    }

    reduce(localRho, minOp<scalar>());
    forceVector_ *= localRho;
}


void Foam::fv::actuatorLineElement::applyForceField
(
    volVectorField& forceField
)
{
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
        Info<< "    sphereRadius: " << sphereRadius << endl;
    }
}


void Foam::fv::actuatorLineElement::calculateInflowVelocity
(
    const volVectorField& Uin
)
{
    // Find local flow velocity by interpolating to element location
    inflowVelocity_ = vector(VGREAT, VGREAT, VGREAT);
    vector inflowVelocityPoint = position_;
    interpolationCellPoint<vector> UInterp(Uin);
    
    // If the flow only is sampled in the center
    if (velocitySampleRadius_ <= 0.0)
    {
        label inflowCellI = findCell(inflowVelocityPoint);
        if (inflowCellI >= 0)
        {
            inflowVelocity_ = UInterp.interpolate
            (
                inflowVelocityPoint,
                inflowCellI
            );
        }

        // Reduce inflow velocity over all processors
        reduce(inflowVelocity_, minOp<vector>());
    }
    // If the flow is sampled by using a circle around position_
    else
    {
        // Circle radius should be normalized with epsilon
        scalar sampleRadius = calcProjectionEpsilon()*velocitySampleRadius_;

        // Unit vector in chordwise direction
        vector chordNormal = chordDirection_ / mag(chordDirection_);
        
        // Calculate mean value over all circle points
        vector velocitySum = vector(0.0, 0.0, 0.0);
        for (label point = 0; point < nVelocitySamples_; point++)
        {
            vector sampleVelocity = vector(VGREAT, VGREAT, VGREAT);

            // distribute the points evenly in terms of angular distance
            scalar pointAngle = Foam::constant::mathematical::pi * 2.0 * point/
                                nVelocitySamples_;
            scalar chordDist = sampleRadius * Foam::cos(pointAngle);
            scalar normalDist = sampleRadius * Foam::sin(pointAngle);
            vector samplePoint = inflowVelocityPoint +
                                 chordDist * chordNormal +
                                 normalDist * planformNormal_;

            // Sample the velocity
            label sampleCellI = findCell(samplePoint);
            if (sampleCellI >= 0)
            {
                sampleVelocity = UInterp.interpolate
                (
                    samplePoint,
                    sampleCellI
                );
            }

            // Reduce inflow velocity over all processors
            reduce(sampleVelocity, minOp<vector>());

            // If inflow velocity is not detected, position is not in the mesh
            if (not (sampleVelocity[0] < VGREAT))
            {
                // Raise fatal error since inflow velocity cannot be detected
                FatalErrorIn("void actuatorLineElement::calculateForce()")
                    << "Inflow velocity point for " << name_
                    << " not found in mesh"
                    << abort(FatalError);
            }

            velocitySum = velocitySum + sampleVelocity;
        }

        // Set inflow Velocity as the mean value
        inflowVelocity_ = 1.0 / nVelocitySamples_ * velocitySum;
    }

    // If inflow velocity is not detected, position is not in the mesh
    if (not (inflowVelocity_[0] < VGREAT))
    {
        // Raise fatal error since inflow velocity cannot be detected
        FatalErrorIn("void actuatorLineElement::calculateForce()")
            << "Inflow velocity point for " << name_
            << " not found in mesh"
            << abort(FatalError);
    }
}


void Foam::fv::actuatorLineElement::createOutputFile()
{
    fileName dir;

    if (Pstream::parRun())
    {
        dir = mesh_.time().path()/"../postProcessing/actuatorLineElements"
            / mesh_.time().timeName();
    }
    else
    {
        dir = mesh_.time().path()/"postProcessing/actuatorLineElements"
            / mesh_.time().timeName();
    }

    if (not isDir(dir))
    {
        mkDir(dir);
    }

    outputFile_ = new OFstream(dir/name_ + ".csv");

    *outputFile_<< "time,root_dist,x,y,z,rel_vel_mag,Re,alpha_deg,"
                << "alpha_geom_deg,cl,cd,fx,fy,fz,end_effect_factor,"
                << "c_ref_t,c_ref_n,f_ref_t,f_ref_n" << endl;
}


void Foam::fv::actuatorLineElement::writePerf()
{
    scalar time = mesh_.time().value();

    scalar pitch;

    dict_.lookup("pitch") >> pitch;    

    // write time,root_dist,x,y,z,rel_vel_mag,Re,alpha_deg,alpha_geom_deg,cl,cd,
    // fx,fy,fz,end_effect_factor,c_ref_t,c_ref_n,f_ref_t,f_ref_n
    *outputFile_<< time << "," << rootDistance_ << "," << position_.x() << ","
                << position_.y() << "," << position_.z() << ","
                << mag(relativeVelocity_) << "," << Re_ << "," << angleOfAttack_
                << "," << angleOfAttackGeom_ << "," << liftCoefficient_ << ","
                << dragCoefficient_ << "," << forceVector_.x() << ","
                << forceVector_.y() << "," << forceVector_.z() << ","
                << endEffectFactor_ << "," << tangentialRefCoefficient() << ","
                << normalRefCoefficient() << "," << tangentialRefForce() << ","
                << normalRefForce() << "," << dict_.lookup("pitch") << endl;
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
    meshBoundBox_(mesh_.points(), false),
    planformNormal_(vector::zero),
    velocity_(vector::zero),
    forceVector_(vector::zero),
    relativeVelocity_(vector::zero),
    relativeVelocityGeom_(vector::zero),
    angleOfAttack_(0.0),
    angleOfAttackGeom_(0.0),
    liftCoefficient_(0.0),
    dragCoefficient_(0.0),
    momentCoefficient_(0.0),
    profileName_(dict.lookup("profileName")),
    profileData_(profileName_, dict.subDict("profileData"), debug),
    dynamicStallActive_(false),
    omega_(0.0),
    chordMount_(0.25),
    flowCurvatureActive_(false),
    flowCurvatureModelName_("none"),
    velocityLE_(vector::zero),
    velocityTE_(vector::zero),
    writePerf_(false),
    rootDistance_(0.0),
    endEffectFactor_(1.0),
    addedMassActive_(dict.lookupOrDefault("addedMass", false)),
    addedMass_(mesh.time(), dict.lookupOrDefault("chordLength", 1.0), debug)
{
    meshBoundBox_.inflate(1e-6);
    read();
    if (writePerf_)
    {
        createOutputFile();
    }
}

// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

Foam::fv::actuatorLineElement::~actuatorLineElement()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::word& Foam::fv::actuatorLineElement::name() const
{
    return name_;
}


const Foam::scalar& Foam::fv::actuatorLineElement::chordLength() const
{
    return chordLength_;
}


const Foam::scalar& Foam::fv::actuatorLineElement::spanLength()
{
    return spanLength_;
}


const Foam::vector& Foam::fv::actuatorLineElement::position()
{
    return position_;
}


const Foam::vector& Foam::fv::actuatorLineElement::velocity()
{
    return velocity_;
}


const Foam::vector& Foam::fv::actuatorLineElement::relativeVelocity()
{
    return relativeVelocity_;
}


const Foam::vector& Foam::fv::actuatorLineElement::relativeVelocityGeom()
{
    return relativeVelocityGeom_;
}


const Foam::scalar& Foam::fv::actuatorLineElement::angleOfAttack()
{
    return angleOfAttack_;
}


const Foam::scalar& Foam::fv::actuatorLineElement::angleOfAttackGeom()
{
    return angleOfAttackGeom_;
}


const Foam::scalar& Foam::fv::actuatorLineElement::liftCoefficient()
{
    return liftCoefficient_;
}


const Foam::scalar& Foam::fv::actuatorLineElement::dragCoefficient()
{
    return dragCoefficient_;
}


const Foam::scalar& Foam::fv::actuatorLineElement::momentCoefficient()
{
    return momentCoefficient_;
}


Foam::scalar Foam::fv::actuatorLineElement::tangentialRefCoefficient()
{
    return profileData_.convertToCRT
    (
        liftCoefficient_,
        dragCoefficient_,
        inflowRefAngle()
    );
}


Foam::scalar Foam::fv::actuatorLineElement::tangentialRefForce()
{
    return 0.5 * chordLength_ * tangentialRefCoefficient()
        * magSqr(relativeVelocity_);
}


Foam::scalar Foam::fv::actuatorLineElement::normalRefCoefficient()
{
    return profileData_.convertToCRN
    (
        liftCoefficient_,
        dragCoefficient_,
        inflowRefAngle()
    );
}


Foam::scalar Foam::fv::actuatorLineElement::normalRefForce()
{
    return 0.5 * chordLength_ * normalRefCoefficient()
        * magSqr(relativeVelocity_);
}


Foam::scalar Foam::fv::actuatorLineElement::inflowRefAngle()
{
    // Calculate inflow velocity angle in degrees (AFTAL Phi)
    scalar inflowVelAngleRad = acos
    (
        (-relativeVelocity_ & chordRefDirection_)
        / (mag(relativeVelocity_) * mag(chordRefDirection_))
    );
    return radToDeg(inflowVelAngleRad);
}


const Foam::scalar& Foam::fv::actuatorLineElement::rootDistance()
{
    return rootDistance_;
}


void Foam::fv::actuatorLineElement::calculateForce
(
    const volVectorField& Uin
)
{
    scalar pi = Foam::constant::mathematical::pi;

    // Calculate vector normal to chord--span plane
    planformNormal_ = -chordDirection_ ^ spanDirection_;
    planformNormal_ /= mag(planformNormal_);

    if (debug)
    {
        Info<< "Calculating force contribution from actuatorLineElement "
            << name_ << endl;
        Info<< "    position: " << position_ << endl;
        Info<< "    chordDirection: " << chordDirection_ << endl;
        Info<< "    spanDirection: " << spanDirection_ << endl;
        Info<< "    elementVelocity: " << velocity_ << endl;
        Info<< "    planformNormal: " << planformNormal_ << endl;
    }

    // Find local flow velocity by interpolating to element location
    calculateInflowVelocity(Uin);

    // Subtract spanwise component of inflow velocity
    vector spanwiseVelocity = spanDirection_
                            * (inflowVelocity_ & spanDirection_)
                            / magSqr(spanDirection_);
    inflowVelocity_ -= spanwiseVelocity;

    // Calculate relative velocity and Reynolds number
    relativeVelocity_ = inflowVelocity_ - velocity_;
    Re_ = mag(relativeVelocity_)*chordLength_/nu_;

    // Calculate angle of attack (radians)
    scalar angleOfAttackRad = asin((planformNormal_ & relativeVelocity_)
                            / (mag(planformNormal_)
                            *  mag(relativeVelocity_)));
    scalar angleOfAttackUncorrected = radToDeg(angleOfAttackRad);
    relativeVelocityGeom_ = freeStreamVelocity_ - velocity_;
    angleOfAttackGeom_ = asin((planformNormal_ & relativeVelocityGeom_)
                       / (mag(planformNormal_)*mag(relativeVelocityGeom_)));
    angleOfAttackGeom_ *= 180.0/pi;

    // Apply flow curvature correction to angle of attack
    if (flowCurvatureActive_)
    {
        correctFlowCurvature(angleOfAttackRad);
    }

    // Calculate angle of attack in degrees
    angleOfAttack_ = radToDeg(angleOfAttackRad);

    // Update Reynolds number of profile data
    profileData_.updateRe(Re_);

    // Lookup lift and drag coefficients
    lookupCoefficients();

    if (debug)
    {
        Info<< "    inflowVelocity: " << inflowVelocity_ << endl;
        Info<< "    relativeVelocity: " << relativeVelocity_ << endl;
        Info<< "    Reynolds number: " << Re_ << endl;
        Info<< "    Geometric angle of attack (degrees): "
            << angleOfAttackGeom_ << endl;
        Info<< "    Angle of attack (uncorrected, degrees): "
            << angleOfAttackUncorrected << endl;
        Info<< "    Angle of attack (corrected, degrees): "
            << angleOfAttack_ << endl;
    }

    // Correct coefficients with dynamic stall model
    if (dynamicStallActive_)
    {
        dynamicStall_->correct
        (
            mag(relativeVelocity_),
            angleOfAttack_,
            liftCoefficient_,
            dragCoefficient_,
            momentCoefficient_
        );
    }

    // Correct for added mass effects
    if (addedMassActive_)
    {
        addedMass_.correct
        (
            liftCoefficient_,
            dragCoefficient_,
            momentCoefficient_,
            degToRad(angleOfAttack_),
            mag(chordDirection_ & relativeVelocity_),
            mag(planformNormal_ & relativeVelocity_)
        );
    }

    // Apply end effect correction factor to lift coefficient
    liftCoefficient_ *= endEffectFactor_;

    // Calculate force per unit density
    scalar area = chordLength_ * spanLength_;
    scalar magSqrU = magSqr(relativeVelocity_);
    scalar lift = 0.5*area*liftCoefficient_*magSqrU;
    scalar drag = 0.5*area*dragCoefficient_*magSqrU;
    vector liftDirection = relativeVelocity_ ^ spanDirection_;
    liftDirection /= mag(liftDirection);
    vector dragDirection = relativeVelocity_/mag(relativeVelocity_);
    forceVector_ = lift*liftDirection + drag*dragDirection;

    if (debug)
    {
        Info<< "    liftDirection: " << liftDirection << endl;
        Info<< "    dragDirection: " << dragDirection << endl;
        Info<< "    force (per unit density): " << forceVector_ << endl;
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
    if (rotateVelocity)
    {
        velocity_ = RM & velocity_;
        chordRefDirection_ = RM & chordRefDirection_;
    }

    if (debug)
    {
        Info<< "Final position: " << position_ << endl;
        Info<< "Final chordDirection: " << chordDirection_ << endl;
        Info<< "Final chordRefDirection: " << chordRefDirection_ << endl;
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

    scalar angleLE = 0.0;
    scalar angleTE = 0.0;
    if (radius > 0.0)
    {
        // Set velocity at leading edge
        scalar radiusLE = sqrt(magSqr(0.25*chordLength_) + magSqr(radius));
        angleLE = atan2(0.25*chordLength_, radius);
        velocityLE_ = velocity_*radiusLE/radius;
        rotateVector(velocityLE_, vector::zero, spanDirection_, angleLE);

        // Set velocity at trailing edge
        scalar radiusTE = sqrt(magSqr(0.75*chordLength_) + magSqr(radius));
        angleTE = atan2(-0.75*chordLength_, radius);
        velocityTE_ = velocity_*radiusTE/radius;
        rotateVector(velocityTE_, vector::zero, spanDirection_, angleTE);
    }

    // Also set omega for flow curvature correction
    setOmega(omega);

    if (debug)
    {
        Info<< "    Radius: " << radius << endl;
        Info<< "    Final velocity: " << velocity_ << endl;
        Info<< "    Leading edge velocity: " << velocityLE_ << endl;
        Info<< "    Trailing edge velocity: " << velocityTE_ << endl;
        Info<< "    Leading edge velocity angle (radians): "
            << angleLE << endl;
        Info<< "    Trailing edge velocity angle (radians): "
            << angleTE << endl;
    }
}


void Foam::fv::actuatorLineElement::scaleVelocity(scalar scale)
{
    velocity_ *= scale;
}


const Foam::vector& Foam::fv::actuatorLineElement::force()
{
    return forceVector_;
}


Foam::vector Foam::fv::actuatorLineElement::moment(vector point)
{
    // Calculate radius vector
    vector radius = position_ - point;
    vector moment = radius ^ forceVector_;
    vector pitchingMoment = 0.5*chordLength_*chordLength_*spanLength_
                          * momentCoefficient_*magSqr(relativeVelocity_)
                          * spanDirection_;
    return moment + pitchingMoment;
}


void Foam::fv::actuatorLineElement::addSup
(
    fvMatrix<vector>& eqn,
    volVectorField& forceField
)
{
    volVectorField forceFieldI
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
            forceField.dimensions(),
            vector::zero
        )
    );

    const volVectorField& Uin(eqn.psi());
    calculateForce(Uin);
    applyForceField(forceFieldI);

    // Add force to total actuator line force
    forceField += forceFieldI;

    // Write performance to file
    if (writePerf_ and Pstream::master())
    {
        writePerf();
    }
}


void Foam::fv::actuatorLineElement::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    volVectorField& forceField
)
{
    volVectorField forceFieldI
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
            forceField.dimensions()/rho.dimensions(),
            vector::zero
        )
    );

    const volVectorField& Uin(eqn.psi());
    calculateForce(Uin);
    applyForceField(forceFieldI);

    // Multiply force vector by local density
    multiplyForceRho(rho);

    // Multiply this element's force field by density field
    forceFieldI *= rho;

    // Add force to total actuator line force
    forceField += forceFieldI;

    // Write performance to file
    if (writePerf_ and Pstream::master())
    {
        writePerf();
    }
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


void Foam::fv::actuatorLineElement::setEndEffectFactor(scalar factor)
{
    endEffectFactor_ = factor;
}


void Foam::fv::actuatorLineElement::setVelocitySampleRadius(scalar radius)
{
    velocitySampleRadius_ = radius;
}


void Foam::fv::actuatorLineElement::setNVelocitySamples(label nSamples)
{
    nVelocitySamples_ = nSamples;
}


// ************************************************************************* //
