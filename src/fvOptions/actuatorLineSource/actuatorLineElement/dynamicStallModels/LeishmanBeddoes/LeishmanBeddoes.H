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
        
        //- Chord length in meters
        scalar c_;
        
        //- Speed of sound (m/s)
        scalar a_;


public:

    //- Runtime type information
    TypeName("LeishmanBeddoes");


    // Constructors

        //- Construct null
        LeishmanBeddoes(const dictionary& dict, const word& modelName);

        //- Construct as copy
        LeishmanBeddoes(const LeishmanBeddoes&);


    // Selectors

        //- Select null constructed
        static autoPtr<LeishmanBeddoes> New
        (
            const dictionary& dict,
            const word& modelName
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

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
