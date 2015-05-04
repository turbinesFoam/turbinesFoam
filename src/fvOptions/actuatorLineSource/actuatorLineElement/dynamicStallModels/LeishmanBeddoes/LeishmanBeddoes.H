/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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
        
        //- Current time value
        scalar time_;
        
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
        
        
    // Protected member functions
    
        // Evaluate the static foil data, i.e., calculate CNAlpha, CN1, 
        // alpha1, S1, S2, and CD0
        void evalStaticData
        (
            List<scalar> alphaDegList,
            List<scalar> clList,
            List<scalar> cdList
        );
    
        // Update previous time step values
        void update();


public:

    //- Runtime type information
    TypeName("LeishmanBeddoes");


    // Constructors

        //- Construct from components
        LeishmanBeddoes
        (
            const dictionary& dict, 
            const word& modelName,
            const scalar startTime
        );

        //- Construct as copy
        LeishmanBeddoes(const LeishmanBeddoes&);


    // Selectors

        //- Select null constructed
        static autoPtr<LeishmanBeddoes> New
        (
            const dictionary& dict,
            const word& modelName,
            const scalar startTime
        );


    //- Destructor
    ~LeishmanBeddoes();


    // Member Functions

        // Access

        // Check

        // Edit

        // Write
        
        // Evaluation
        
            // Correct lift and drag coefficients
            virtual void correct(scalar alphaDeg, scalar& cl, scalar& cd);
            
            // Correct lift and drag coefficients
            virtual void correct
            (
                scalar time,
                scalar magU,
                scalar alphaDeg,
                scalar& cl,
                scalar& cd,
                List<scalar> alphaDegList,
                List<scalar> clList,
                List<scalar> cdList
            );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
