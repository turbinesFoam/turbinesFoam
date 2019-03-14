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
    Foam::fv::axialFlowTurbineALSource

Description
    Cell based momentum source that represents a axial-flow turbine

SourceFiles
    axialFlowTurbineALSource.C

\*---------------------------------------------------------------------------*/

#ifndef axialFlowTurbineALSource_H
#define axialFlowTurbineALSource_H

#include "turbineALSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                      Class axialFlowTurbineALSource Declaration
\*---------------------------------------------------------------------------*/

class axialFlowTurbineALSource
:
    public turbineALSource
{

protected:

    // Protected data

        //- Switch for hub present
        bool hasHub_;

        //- Switch for tower present
        bool hasTower_;

        //- Switch for nacelle present
        bool hasNacelle_;

        //- Actuator line for hub
        autoPtr<actuatorLineSource> hub_;

        //- Actuator line for tower
        autoPtr<actuatorLineSource> tower_;

        //- Actuator line for nacelle
        autoPtr<actuatorLineSource> nacelle_;

        //- Hub dictionary
        dictionary hubDict_;

        //- Tower definitions dictionary
        dictionary towerDict_;

        //- Nacelle dictionary
        dictionary nacelleDict_;

        //- Vertical direction
        vector verticalDirection_;

        //- Initial azimuthal direction
        vector azimuthalDirection_;

        //- Switch to include tower drag in turbine drag coefficient
        bool includeTowerDrag_;

        //- Switch to include nacelle drag in turbine drag coefficient
        bool includeNacelleDrag_;

        //- End effects dictionary
        dictionary endEffectsDict_;

        //- Switch to activate rotor-level end effects model
        bool endEffectsActive_;

        //- End effects model name
        word endEffectsModel_;


    // Protected Member Functions

        //- Create the co-ordinate system
        virtual void createCoordinateSystem();

        //- Create actuator lines for blades
        virtual void createBlades();

        //- Create actuator line for hub
        void createHub();

        //- Create actuator line for tower
        void createTower();

        //- Create actuator line for nacelle
        void createNacelle();

        //- Calculate end end effects at rotor level
        void calcEndEffects();

        //- Use turbine base class rotate method
        using turbineALSource::rotate;

        //- Rotate the turbine a specified angle about its axis
        virtual void rotate(scalar radians);

        //- Yaw the turbine, i.e., rotate about vertical axis
        void yaw(scalar radians);


public:

    //- Runtime type information
    TypeName("axialFlowTurbineALSource");


    // Constructors

        //- Construct from components
        axialFlowTurbineALSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~axialFlowTurbineALSource();


    // Member Functions

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
