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
    Foam::fv::crossFlowTurbineALSource

Description
    Cell based momentum source that represents a cross-flow turbine

SourceFiles
    crossFlowTurbineALSource.C

\*---------------------------------------------------------------------------*/

#ifndef crossFlowTurbineALSource_H
#define crossFlowTurbineALSource_H

#include "turbineALSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                      Class crossFlowTurbineALSource Declaration
\*---------------------------------------------------------------------------*/

class crossFlowTurbineALSource
:
    public turbineALSource
{

protected:

    // Protected data

        //- Switch for struts
        bool hasStruts_;

        //- Switch for shaft
        bool hasShaft_;

        //- List of actuator lines as struts
        PtrList<actuatorLineSource> struts_;

        //- Actuator line for shaft
        autoPtr<actuatorLineSource> shaft_;

        //- Strut definitions dictionary
        dictionary strutsDict_;

        //- Shaft dictionary
        dictionary shaftDict_;


    // Protected Member Functions

        //- Create actuator lines for blades
        virtual void createBlades();

        //- Create actuator lines for struts
        void createStruts();

        //- Create actuator line for shaft
        void createShaft();

        //- Create the co-ordinate system
        virtual void createCoordinateSystem();

        //- Use turbine base class rotate method
        using turbineALSource::rotate;

        //- Rotate the turbine a specified angle about its axis
        virtual void rotate(scalar radians);


public:

    //- Runtime type information
    TypeName("crossFlowTurbineALSource");


    // Constructors

        //- Construct from components
        crossFlowTurbineALSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~crossFlowTurbineALSource();


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
