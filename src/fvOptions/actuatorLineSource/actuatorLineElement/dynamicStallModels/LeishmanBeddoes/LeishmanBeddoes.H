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
    Foam::fv::LeishmanBeddoes

Description
    Leishman-Beddoes dynamic stall model for use with actuatorLineElement.

SourceFiles
    LeishmanBeddoes.C

\*---------------------------------------------------------------------------*/

#ifndef LeishmanBeddoes_H
#define LeishmanBeddoes_H

#include "dynamicStallModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                         Class LeishmanBeddoes Declaration
\*---------------------------------------------------------------------------*/

class LeishmanBeddoes
:
    public dynamicStallModel
{

protected:

    // Protected data

        //- Chord length in meters
        scalar c_;

        //- Normal force coefficient from circulatory effects
        scalar CNC_;

        //- Equivalent angle of attack in radians
        scalar alphaEquiv_;

        //- Deficiency function for equivalent angle of attack
        scalar X_;

        //- Previous value of X
        scalar XPrev_;

        //- Deficiency function for equivalent angle of attack
        scalar Y_;

        //- Previous value of Y
        scalar YPrev_;

        //- Mach number
        scalar M_;

        //- Constant in angle of attack deficiency function; 0.3 by default
        scalar A1_;

        //- Constant in angle of attack deficiency function; 0.7 by default
        scalar A2_;

        //- Constant in angle of attack deficiency function; 0.14 by default
        scalar b1_;

        //- Constant in angle of attack deficiency function; 0.53 by default
        scalar b2_;

        //- Angle of attack (radians)
        scalar alpha_;

        //- Previous angle of attack (radians)
        scalar alphaPrev_;

        //- Change in angle of attack (radians)
        scalar deltaAlpha_;

        //- Previous value of deltaAlpha
        scalar deltaAlphaPrev_;

        //- Time step in seconds
        scalar deltaT_;

        //- Nondimensional time step
        scalar deltaS_;

        //- Normal force coefficient from impulsive effects
        scalar CNI_;

        //- Time constant for impulsive loading = c/a
        scalar TI_;

        //- Speed of sound (m/s)
        scalar a_;

        //- Previous time value
        scalar timePrev_;

        //- Deficiency function for the impulsive force response
        scalar D_;

        //- Previous value of D
        scalar DPrev_;

        //- Deficiency function for delayed pressure response
        scalar DP_;

        //- Previous value of DP
        scalar DPPrev_;

        //- Total normal force coefficient under the unsteady attached flow
        scalar CNP_;

        //- Previous value of CNP
        scalar CNPPrev_;

        //- Lagged normal force coefficient
        scalar CNPrime_;

        //- Normal force coefficient slope
        scalar CNAlpha_;

        //- Critical normal force coefficient
        scalar CN1_;

        //- Delayed separation point
        scalar fPrime_;

        //- Previous value of fPrime
        scalar fPrimePrev_;

        //- Second order delayed separation point
        scalar fDoublePrime_;

        //- Delayed angle of attack (radians)
        scalar alphaPrime_;

        //- Angle of attack for the breakpoint of separation (radians)
        scalar alpha1_;

        //- Deficiency function for the delayed separation point
        scalar DF_;

        //- Previous value of DF
        scalar DFPrev_;

        //- Normal force coefficient for the delayed separation point
        scalar CNF_;

        //- Vortex lift
        scalar CV_;

        //- Previous value of CV
        scalar CVPrev_;

        //- Vortex lift contribution to normal force coefficient
        scalar CNV_;

        //- Previous value of CNV
        scalar CNVPrev_;

        //- Total normal force coefficient
        scalar CN_;

        //- Tangential force coefficient
        scalar CT_;

        //- Tangential force efficiency factor
        scalar eta_;

        //- Separation point calculation constants
        scalar S1_;
        scalar S2_;

        //- Drag coefficient at zero angle of attack
        scalar CD0_;

        //- Switch to indicate stall
        bool stalled_;

        //- Previous stall state
        bool stalledPrev_;

        //- Time constant associated with the pressure response, 1.7 by default
        scalar Tp_;

        //- Time constant associated with the separation point, 3 by default
        scalar Tf_;

        //- Time constant associated with the vortex lift
        scalar Tv_;

        //- Time constant associated with the vortex convection
        scalar Tvl_;

        //- Nondimensional vortex tracking time
        scalar tau_;

        //- Previous value of tau
        scalar tauPrev_;

        //- Number of new times evaluated -- does not increment for subcycles
        label nNewTimes_;

        //- Static stall angle (radians)
        scalar alphaSS_;

        //- Critical trailing edge separation point
        scalar fCrit_;

        //- Magnitude of the relative velocity
        scalar magU_;

        //- Offset of aerodynamic center (0.25 - x_ac)
        scalar K0_;

        //- Linear pitching moment curve fit coefficient
        scalar K1_;

        //- Sinusoidal pitching moment curve fit coefficient
        scalar K2_;

        //- Exponent on f in moment coefficient fit equation
        scalar cmFitExponent_;

        //- Total moment coefficient
        scalar CM_;

        //- Reynolds number
        scalar Re_;


    // Protected member functions

        //- Calculate the equivalent angle of attack
        virtual void calcAlphaEquiv();

        //- Evaluate the static foil data, i.e., calculate CNAlpha, CN1,
        // alpha1, and CD0
        virtual void evalStaticData();

        //- Interpolate static profile data based on Reynolds number
        scalar interpolateStaticParam(word paramName);

        //- Calculate unsteady quantities
        virtual void calcUnsteady();

        //- Calculate S1 and S2 constants
        virtual void calcS1S2
        (
            scalar B=0.3,
            scalar C=0.04,
            scalar D=0.66
        );

        //- Calculate K1 and K2 coefficients
        virtual void calcK1K2();

        //- Calculate separated flow quantities
        virtual void calcSeparated();

        //- Convert normal force coefficient to trailing edge separation point
        virtual List<scalar> cnToF
        (
            List<scalar> cnList,
            List<scalar> alphaRadList,
            bool limit=true
        );

        //- Update previous time step values
        virtual void update();


public:

    //- Runtime type information
    TypeName("LeishmanBeddoes");


    // Constructors

        //- Construct from components
        LeishmanBeddoes
        (
            const dictionary& dict,
            const word& modelName,
            const Time& time,
            profileData& profileData
        );


    // Selectors

        //- Select from components
        static autoPtr<LeishmanBeddoes> New
        (
            const dictionary& dict,
            const word& modelName,
            const Time& time,
            profileData& profileData
        );


    //- Destructor
    ~LeishmanBeddoes();


    // Member Functions

        //- Correct lift, drag, and moment coefficients
        virtual void correct
        (
            scalar magU,
            scalar alphaDeg,
            scalar& cl,
            scalar& cd,
            scalar& cm
        );


        // Parallel running

            //- Reduce to set data equal on all processors
            virtual void reduceParallel(bool inMesh);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
