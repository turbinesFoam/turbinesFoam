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
    Foam::fv::LeishmanBeddoesSGC

Description
    Leishman-Beddoes dynamic stall model modified by Sheng et al. (2008).

SourceFiles
    LeishmanBeddoesSGC.C

\*---------------------------------------------------------------------------*/

#ifndef LeishmanBeddoesSGC_H
#define LeishmanBeddoesSGC_H

#include "LeishmanBeddoes3G.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                         Class LeishmanBeddoesSGC Declaration
\*---------------------------------------------------------------------------*/

class LeishmanBeddoesSGC
:
    public LeishmanBeddoes3G
{

protected:

    // Protected data

        //- Time constant associated with delayed angle of attack
        scalar TAlpha_;

        //- Previous value of alphaPrime
        scalar alphaPrimePrev_;

        //- Critical angle of attack
        scalar alphaCrit_;

        //- Dynamic stall angle of attack for high pitch rate
        scalar alphaDS0_;

        //- Reduced pitch rate
        scalar r_;

        //- Threshold value of reduced pitch rate
        scalar r0_;

        //- Constant associated with vortex lift
        scalar B1_;

        //- Constant associated with vortex moment
        scalar B2_;

        //- Constant associated with the tangential force
        scalar E0_;

        //- Difference between alphaDS0 and alphaSS (degrees)
        scalar alphaDS0DiffDeg_;

        //- Deficiency function for alphaPrime
        scalar DAlpha_;

        //- Previous value of DAlpha
        scalar DAlphaPrev_;


    // Protected member functions

        //- Calculate unsteady quantities
        virtual void calcUnsteady();

        //- Calculate separated flow quantities
        virtual void calcSeparated();

        //- Update previous time step values
        virtual void update();


public:

    //- Runtime type information
    TypeName("LeishmanBeddoesSGC");


    // Constructors

        //- Construct from components
        LeishmanBeddoesSGC
        (
            const dictionary& dict,
            const word& modelName,
            const Time& time,
            profileData& profileData
        );


    // Selectors

        //- Select from components
        static autoPtr<LeishmanBeddoesSGC> New
        (
            const dictionary& dict,
            const word& modelName,
            const Time& time,
            profileData& profileData
        );


    //- Destructor
    ~LeishmanBeddoesSGC();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
