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

        
    // Protected member functions
    
        // Calculate the equivalent angle of attack
        virtual void calcAlphaEquiv();
        
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
    TypeName("LeishmanBeddoesSGC");


    // Constructors

        //- Construct from components
        LeishmanBeddoesSGC
        (
            const dictionary& dict, 
            const word& modelName,
            const Time& time
        );


    // Selectors

        //- Select null constructed
        static autoPtr<LeishmanBeddoesSGC> New
        (
            const dictionary& dict,
            const word& modelName,
            const Time& time
        );


    //- Destructor
    ~LeishmanBeddoesSGC();


    // Member Functions

        // Access

        // Check

        // Edit

        // Write
        
        // Evaluation

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
