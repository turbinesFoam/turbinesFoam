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

Class
    Foam::fv::actuatorLineSource

Description
    Actuator line class, which is a collection of actuator line elements.

SourceFiles
    actuatorLineSource.C

\*---------------------------------------------------------------------------*/

#ifndef actuatorLineSource_H
#define actuatorLineSource_H

#include "List.H"
#include "dictionary.H"
#include "vector.H"
#include "actuatorLineElement.H"
#include "cellSetOption.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                         Class actuatorLineSource Declaration
\*---------------------------------------------------------------------------*/

class actuatorLineSource
:
    public cellSetOption
{

protected:

    // Protected data

        //- Total length of geometry in meters
        scalar totalLength_;

        //- Dictionary for the profiles' sectional coefficient data
        dictionary profileData_;

        //- Element geometry
        // point, spanDirection, chordLength, chordRefDir, pitch
        List<List<List<scalar> > > elementGeometry_;

        //- Element profiles
        List<word> elementProfiles_;

        //- Number of blade elements
        label nElements_;

        //- question //Velocity of the element associated to prescribed/rigidBody motion
        // vector motionVelocity_;

        //- Free stream velocity
        vector freeStreamVelocity_;

        //- Free stream direction
        vector freeStreamDirection_;

        //- Pitch relative flow velocity
        vector pitchRelativeVelocity_;

        //- Total force vector from all elements
        vector force_;

        //- Force field from all elements
        volVectorField forceField_;

        //- List of actuator line elements
        PtrList<actuatorLineElement> elements_;

        //- Switch for writing performance
        bool writePerf_;

        //- Output file stream
        OFstream* outputFile_;

        //- Switch for harmonic pitching
        bool harmonicPitchingActive_;

        // New- Switch for wind pitching
        bool windPitchingActive_;

        // New
        scalar lastPitchAngle_;

        //- Reduced frequency of harmonic pitching = omega*c/(2*UInfty)
        scalar reducedFreq_;

        //- Amplitude of harmonic pitching in degrees
        scalar pitchAmplitude_;

        // - Height of turbine's hub
        scalar hubHeight_;

        //- Time value to track whether to move
        scalar lastMotionTime_;

        //- Mean chord length of all elements
        scalar chordLength_;

        //- Relative velocity at the hub
        scalar Urel_;

        //- Aspect ratio
        scalar aspectRatio_;

        //- Switch for correcting end effects
        bool endEffectsActive_;

        //- Switch for multiphase simulations
        bool multiPhase_;

        //- Name of the phase which properties are to be used
        word phaseName_;

        //- Prescribed motion dictionary
        dictionary prescribedMotionDict_;

        //- Switch for prescribed motion
        bool prescribedMotionActive_;

        //- Type of prescribed motion
        word prescribedMotionType_;

        //- Prescribed motion coefficients dictionary
        dictionary prescribedMotionCoeffs_;

        //- Rigid body motion dictionary
        dictionary rigidBodyMotionDict_;

        //- Switch for rigid body motion
        bool rigidBodyMotionActive_;

        //- Dictionary from "sixDoFRigidBodyMotionState"
        dictionary rigidBodyDict_;

        //- Orientation at t=0 of the prescribed/rigidBody motion
        vector rot0_;

        //- Center of rotation of the prescribed/rigidBody motion
        vector rotCenter_;

        //- rigidBody linaer velocity
        vector linearVelocity_;

        //- rigidBody angular velocity
        vector angularVelocity_;

        //- Center of rotation at previous time of the prescribed/rigidBody motion
        vector prevRotCenter_;

        //- Floater's orientation as a zy'x'' rotation matrix
        // from floater (F) to inertial (I) frame
        tensor orientation_;

        //- Floater's orientation at previous time as a zy'x'' rotation Matrix
        tensor prevOrientation_;

        //- Orientation from rigid body library
        tensor rBOrientation_;

        //- For rigid body motion, orientation aligned with inertia principal axes
        tensor refOrientation_;

        //- False if the AL is defined wrt to orientation_ = I, and thus must be 
        // rotated to match rigid body orientation before the first timestep.
        // True would mean that the user has defined the AL aligned with
        // the rigid body initial orientation.  
        bool rigidBodyAligned_;

    // Protected Member Functions

        //- Create actuator line elements
        void createElements();

        //- Read dictionary
        bool read(const dictionary& dict);

        //- Create the performance output file
        virtual void createOutputFile();

        //- Write performance to CSV
        void writePerf();

        //- Calculate end effects from lifting line theory
        void calcEndEffects();

        //- Execute harmonic pitching for a single time step
        void harmonicPitching();

        //-New Execute wind pitching for a single time step
        void windPitching(const volVectorField& Uin);

        //- Initialise motion orientation according to prescribed motion
        void prescribedMotionInitialise();

        //- Compute specific motion parameters according to motionType
        void prescribedMotionParameters(vector &translation, vector &rotation, vector &velocity, vector &omega, const scalar t, const scalar dt);

        //- Read rigid body dictionary
        void readRigidBodyDict(const fvMesh& mesh);

        //- Initialise rigid body
        void rigidBodyInitialise(const fvMesh& mesh);

public:

    //- Runtime type information
    TypeName("actuatorLineSource");

    // Selectors

    //- Return a reference to the selected fvOption model
    static autoPtr<actuatorLineSource> New
    (
        const word& name,
        const dictionary& dict,
        const fvMesh& mesh
    );


    //- Constructor
    actuatorLineSource
    (
        const word& name,
        const word& modelType,
        const dictionary& dict,
        const fvMesh& mesh
    );


    //- Destructor
    virtual ~actuatorLineSource();


    // Member functions

        // Access

            //- Return const reference to the total force vector
            const vector& force();

            //- Return const reference to the force field
            const volVectorField& forceField();

            //- Return reference to element pointer list
            PtrList<actuatorLineElement>& elements();

            //- Return const reference to the current orientation
            const tensor& orientation();

            //- Return pitch relative velocity
            const vector& pitchRelativeVelocity();

            //- Return const reference to the previous orientation
            const tensor& prevOrientation();

            //- Return translation from the prescribed/rigidBody motion at current timestep
            const vector motionTranslation();

            //- Return const reference to the current rotation center
            const vector& rotCenter();

            //- Return const reference to the previous rotation center
            const vector& prevRotCenter();

        // Edit

            //- Translate the actuator line
            void translate(vector translation);

            //- Rotate the actuator line about a specified axis
            void rotate(vector rotationPoint, vector axis, scalar radians);

            //- Rotate the actuator line using a rotation matrix
            void rotate(const vector &rotationPoint, const tensor &rotMatrix);

            //- Pitch the blade about its chord mount
            void pitch(scalar radians);

            //- Pitch the blade about some fraction of the chord
            void pitch(scalar radians, scalar chordFraction);

            //- Move the actuator line according to the prescribed/rigidBody motion. Account for translation velocity.
            void applyMotion(const vector &translation, const vector &velocity, const vector &omega);

            //- Return the intrinsic rotation matrix (zy'x'') given the roll-pitch-yaw angles
            tensor rotAngles2Matrix(const vector &rotAngles);

            //- Tranform a given axis-angle rotation into rotation matrix
            tensor rotAxis2Matrix(const vector &axis, const scalar &angle);

            //- Tranform a given rotation matrix into an axis-angle rotation
            void rotMatrix2Axis(const tensor &rotMatrix, vector &axis, scalar &angle);

            //- Scale the velocity of all elements
            void scaleVelocity(scalar scale);

            //- Set speed of all elements based on rotation
            void setSpeed(vector point, vector axis, scalar omega);

            //- Set prescribedMotion/rigidBody velocity to all elements
            void setMotionVelocity(const vector &velocity);

            //- Add velocity to all elements based on prescribedMotion/rigidBody translation
            void addMotionVelocity(const vector &velocity);

            //- Add velocity to all elements based on prescribedMotion/rigidBody rotation
            void addMotionOmega(const vector &point, const vector &axis, const scalar &omega);

            //- Set rotational velocity of all elements for flow curvature
            //  correction
            void setOmega(scalar omega);

             //- Execute prescribed motion for a single time step
            void prescribeMotion();

            //- Execute rigid body motion
            void rigidBodyMotion(const fvMesh& mesh);

            //- Rotate wrt the rotation center from prescribedMotion/rigidBody (with the orientation matrix)
            void motionRotate();

            //- Rotate according to prescribedMotion/rigidBody initial orientation
            void motionAlign();
            
        // Evaluation

            //- Compute the moment about a given point
            vector moment(vector point);


        // IO

            //- Print dictionary values
            virtual void printCoeffs() const;


        // Source term addition

            //- Source term to momentum equation
            virtual void addSup
            (
                fvMatrix<vector>& eqn,
                const label fieldI
            );

            //- Source term to turbulence scalars
            virtual void addSup
            (
                fvMatrix<scalar>& eqn,
                const label fieldI
            );

            //- Source term to compressible momentum equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                const label fieldI
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
