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
    
        //- Time constant for leading edge pressure response to change in C_l
        scalar timeConstP_;
        
        //- Time constant for trailing edge separation point travel
        scalar timeConstSepPointTravel_;
        
        //- Time constant for leading edge vortex lift indicial function
        scalar timeConstLEVortexLift_;
        
        //- Characteristic leading edge vortex travel time
        scalar timeConstLEVortexTravel_;
        
        //- TODO: Define this
        scalar dp_;
        
        //- TODO: Define this
        scalar dF_;

        //- TODO: Define this
        scalar dCNv_;
        
        //- Switch for leading edge separation
        bool leadingEdgeSeparated_;
        
        //- TODO: Define this
        scalar sLEv_;
        
        //- Previous reference lift coefficient
        scalar liftCoeffRefPrev_;
        
        //- Previous leading edge reference lift coefficient
        scalar liftCoeffRefLEPrev_;
        
        //- TODO: Define this
        scalar FStatePrev_;
        
        //- TODO: Define this
        scalar cvPrev_;

        //- Time step normalized as ds = 2*U*dt/c where U is the relative 
        //  velocity and c is the chord (note: at the same t, s is different 
        //  on each element)
        scalar ds_;

        //- Reference lift coefficient
        scalar liftCoeffRef_;
        
        //- Reference leading edge lift coefficient
        scalar liftCoeffRefLE_;

        //- TODO: Define this
        scalar liftCoeffCritP_;
        
        //- TODO: Define this
        scalar liftCoeffCritN_;

        //- Lift coefficient rate switch
        bool liftCoeffRateFlag_;

        //- TODO: Define this
        scalar FState_;

        //- TODO: Define this
        scalar F_;

        //- TODO: Define this
        scalar cv_;

        //- TODO: Define this
        scalar dcv_;
        
        //- TODO: Define this
        scalar angleOfAttack0_;
        
        //- TODO: Define this
        scalar liftCoeffIdeal_;


public:


    // Constructors

        //- Construct null
        LeishmanBeddoes();

        //- Construct as copy
        LeishmanBeddoes(const LeishmanBeddoes&);


    // Selectors

        //- Select null constructed
        static autoPtr<LeishmanBeddoes> New();


    //- Destructor
    ~LeishmanBeddoes();


    // Member Functions

        // Access

        // Check

        // Edit

        // Write

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
