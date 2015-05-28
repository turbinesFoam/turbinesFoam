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
    Foam::fv::LeishmanBeddoes3G

Description
    Leishman-Beddoes 3rd generation dynamic stall model for use with 
    actuatorLineElement.

SourceFiles
    LeishmanBeddoes3G.C

\*---------------------------------------------------------------------------*/

#ifndef LeishmanBeddoes3G_H
#define LeishmanBeddoes3G_H

#include "LeishmanBeddoes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                         Class LeishmanBeddoes3G Declaration
\*---------------------------------------------------------------------------*/

class LeishmanBeddoes3G
:
    public LeishmanBeddoes
{

protected:

    // Protected data
    
        
    // Protected member functions
    
        // Calculate the equivalent angle of attack
        virtual void calcAlphaEquiv();
    
        // Evaluate the static foil data, i.e., calculate CNAlpha, CN1, 
        // alpha1, and CD0
        virtual void evalStaticData
        (
            List<scalar> alphaDegList,
            List<scalar> clList,
            List<scalar> cdList
        );
        
        // Calculate unsteady quantities
        virtual void calcUnsteady();
        
        // Calculate S1 and S2 constants
        virtual void calcS1S2
        (
            List<scalar> alphaDegList,
            List<scalar> clList,
            List<scalar> cdList
        );
        
        // Calculate separated flow quantities
        virtual void calcSeparated();
    
        // Update previous time step values
        virtual void update();


public:

    //- Runtime type information
    TypeName("LeishmanBeddoes3G");


    // Constructors

        //- Construct from components
        LeishmanBeddoes3G
        (
            const dictionary& dict, 
            const word& modelName,
            const Time& time
        );


    // Selectors

        //- Select null constructed
        static autoPtr<LeishmanBeddoes3G> New
        (
            const dictionary& dict,
            const word& modelName,
            const Time& time
        );


    //- Destructor
    ~LeishmanBeddoes3G();


    // Member Functions

        // Access

        // Check

        // Edit

        // Write
        
        // Evaluation
            
            // Correct lift and drag coefficients
            virtual void correct
            (
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
