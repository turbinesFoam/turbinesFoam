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
    Foam::fv::turbineALSource

Description
    Cell based momentum source, which is a collection of actuatorLineSources

SourceFiles
    turbineALSource.C

\*---------------------------------------------------------------------------*/

#ifndef turbineALSource_H
#define turbineALSource_H

#include "cellSetOption.H"
#include "NamedEnum.H"
#include "actuatorLineSource.H"
#include "volFieldsFwd.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                      Class turbineALSource Declaration
\*---------------------------------------------------------------------------*/

class turbineALSource
:
    public cellSetOption
{

protected:

    // Protected data

        //- Runtime pointer
        const Time& time_;

        //- Time value used to determine whether to rotate
        scalar lastRotationTime_;

        //- Origin of turbine coordinate system
        vector origin_;

        //- Turbine axis of rotation
        vector axis_;

        //- Reference density for incompressible case
        scalar rhoRef_;

        //- Rotational speed in rad/s
        //  Positive anti-clockwise when looking along axis direction
        scalar omega_;

        //- Azimuthal angle in degrees (does not reset each rotation)
        scalar angleDeg_;

        //- Number of blades
        label nBlades_;

        //- Inlet velocity for specified inflow
        vector freeStreamVelocity_;

        //- Free stream velocity direction
        vector freeStreamDirection_;

        //- (Initial) radial direction
        vector radialDirection_;

        //- Rotor tip speed ratio
        scalar tipSpeedRatio_;

        //- Rotor radius for calculating tip speed ratio
        scalar rotorRadius_;

        //- List of actuator lines as blades
        PtrList<actuatorLineSource> blades_;

        //- Blade definitions dictionary
        dictionary bladesDict_;

        //- List of blade names
        wordList bladeNames_;

        //- Dictionary for profile coefficient data
        dictionary profileData_;

        //- Force vector
        vector force_;

        //- Force field (per unit density)
        volVectorField forceField_;

        //- Torque about the axis
        scalar torque_;

        //- Frontal area in m^2
        scalar frontalArea_;

        //- Nondimensional power coefficient
        scalar powerCoefficient_;

        //- Nondimensional drag coefficient
        scalar dragCoefficient_;

        //- Nondimensional torque coefficient
        scalar torqueCoefficient_;

        //- Output file stream
        OFstream* outputFile_;

        //- Dynamic stall dictionary
        dictionary dynamicStallDict_;

        //- Amplitude of prescribed TSR oscillation
        scalar tsrAmplitude_;

        //- Angle of first peak in TSR oscillation (rad)
        scalar tsrPhase_;

        //- Mean tip speed ratio
        scalar meanTSR_;

        //- Individual blade moments about turbine origin
        List<vector> bladeMoments_;


    // Protected Member Functions

        //- Rotate a vector
        void rotateVector
        (
            vector& vectorToRotate,
            vector rotationPoint,
            vector axis,
            scalar radians
        );

        //- Create actuator lines for blades
        virtual void createBlades();

        //- Create the co-ordinate system
        virtual void createCoordinateSystem();

        //- Update tip speed ratio and angular velocity
        virtual void updateTSROmega();

        //- Rotate the turbine according to its angular velocity and time step
        virtual void rotate();

        //- Rotate the turbine a specified angle about its axis
        virtual void rotate(scalar radians);

        //- Print performance
        virtual void printPerf();


public:

    //- Runtime type information
    TypeName("turbineALSource");


    // Constructors

        //- Construct from components
        turbineALSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~turbineALSource();


    // Member Functions

        // Access

            //- Return const access to runTime
            inline const Time& time() const;

            //- Return the reference density for incompressible case
            inline scalar rhoRef() const;


        // Source term addition

            //- Add source term to momentum equation
            virtual void addSup
            (
                fvMatrix<vector>& eqn,
                const label fieldI
            );

            //- Add source term to compressible momentum equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                const label fieldI
            );

            //- Add source term to turbulence model equation
            virtual void addSup
            (
                fvMatrix<scalar>& eqn,
                const label fieldI
            );


        // I-O

            //- Create the turbine output file
            virtual void createOutputFile();

            //- Write the turbine performance to file
            virtual void writePerf();

            //- Write the source properties
            virtual void writeData(Ostream&) const;

            //- Read source dictionary
            virtual bool read(const dictionary& dict);

            //- Print dictionary values
            virtual void printCoeffs() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
