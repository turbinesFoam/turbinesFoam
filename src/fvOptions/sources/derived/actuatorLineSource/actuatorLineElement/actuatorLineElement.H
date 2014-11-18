/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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
    Foam::fv::actuatorLineElement

Description
    An object, e.g. an airfoil, that imparts a force on the flow field 
    calculated via lift and drag coefficients

SourceFiles
    actuatorLineElement.C

\*---------------------------------------------------------------------------*/

#ifndef actuatorLineElement_H
#define actuatorLineElement_H

#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "dictionary.H"
#include "vector.H"
#include "localAxesRotation.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                        Class actuatorLineElement Declaration
\*---------------------------------------------------------------------------*/

class actuatorLineElement
{

protected:

    // Protected data

        //- Coefficients dictionary
        const dictionary dict_;

        //- Name of profile model
        const word name_;

        //- File name (optional)
        fileName fName_;
        
        //- Chord direction
        vector chord_;
        
        //- Chord length
        scalar chordLength_;
        
        //- Element span, for calculating force
        scalar span_;
        
        //- Location of element's half chord
        vector position_;
        
        //- Velocity of the element
        vector velocity_;
        
        //- Fluid velocity near the element
        vector fluidVelocity_;
        
        //- Fluid force on element
        vector force_;
        
        //- Reference density for incompressible case
        scalar rhoRef_;
        
        //- Angle of attack
        scalar angleOfAttack_;
        
        //- Lift coefficient
        scalar liftCoefficient_;
        
        //- Drag coefficient;
        scalar dragCoefficient_;


    // Protected Member Functions

        //- Return true if file name is set
        bool readFromFile() const;


public:

    //- Runtime type information
    TypeName("actuatorLineElement");


        // Declare run-time constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            actuatorLineElement,
            dictionary,
            (
                const dictionary& dict,
                const word& modelName
            ),
            (dict, modelName)
        );


    // Selectors

        //- Return a reference to the selected fvOption model
        static autoPtr<actuatorLineElement> New(const dictionary& dict);


    //- Constructor
    actuatorLineElement
	(
		const word& name,
		const dictionary& dict
	);


    //- Destructor
    virtual ~actuatorLineElement();


    // Member functions

        // Access

            //- Return const access to the element name
            const word& name() const;
            
            //- Return const access to the element location
            const vector& position() const;
            
            //- Return force vector from element
            vector& force();
            
        // Manipulation
        
            //- Translate the element's location
            void translate(vector translationVector);
            
            //- Rotate the element about some axis
            void rotate(vector axis, scalar radians);
            
            //- Pitch the element about its half chord
            void pitch(scalar radians);

        // Evaluation
            
            //- Calculate forces
            void calculate();
            
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
