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
#include "fvMesh.H"
#include "fvMatrices.H"
#include "dynamicStallModel.H"

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
        
        //- Reference to the mesh
        const fvMesh& mesh_;
        
        //- Chord direction
        vector chordDirection_;
        
        //- Chord length
        scalar chordLength_;
        
        //- Span direction
        vector spanDirection_;
        
        //- Element span length, for calculating force
        scalar spanLength_;
        
        //- Location of element's half chord
        vector position_;
        
        //- Velocity of the element
        vector velocity_;
        
        //- Fluid velocity near the element
        vector fluidVelocity_;
        
        //- Free stream direction
        vector freeStreamDirection_;
        
        //- Fluid force on element; a vector that is perpendicular to 
        // spanDirection, 
        vector forceVector_;
        
        //- Reference density for incompressible case
        scalar rhoRef_;
        
        //- Relative flow velocity
        vector relativeVelocity_;
        
        //- Angle of attack (degrees)
        scalar angleOfAttack_;
        
        //- Lift coefficient
        scalar liftCoefficient_;
        
        //- Drag coefficient;
        scalar dragCoefficient_;
        
        //- List of triplets (AoA, C_l, C_d)
        List<List<scalar> > coefficientData_;
        
        //- List of angles of attack
        List<scalar> angleOfAttackList_;
        
        //- List of lift coefficients
        List<scalar> liftCoefficientList_;
        
        //- List of drag coefficients
        List<scalar> dragCoefficientList_;
        
        //- Dynamic stall model
        autoPtr<dynamicStallModel> dynamicStall_;
        
        //- Switch for existence of dynamic stall model
        bool dynamicStallActive_;
        
        //- Kinematic viscosity (for calculating Re)
        scalar nu_;
        
        //- Chord Reynolds number based on relative velocity
        scalar Re_;
        
        //- Angular velocity for flow curvature correction
        scalar omega_;
        
        //- Nondimensional chord mount distance for flow curvature correction
        scalar chordMount_;


    // Protected Member Functions
    
        //- Interpolate a scalar value
        scalar interpolate
        (
            scalar xNew, 
            List<scalar>& xOld, 
            List<scalar>& yOld
        );
        
        //- Lookup force coefficients
        void lookupCoefficients();
        
        //- Calculate projection width epsilon
        scalar calcProjectionEpsilon();


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

        //- Return a reference
        static autoPtr<actuatorLineElement> New(const dictionary& dict);


    //- Constructor
    actuatorLineElement
	(
		const word& name,
		const dictionary& dict,
        const fvMesh& mesh
	);


    //- Destructor
    virtual ~actuatorLineElement();


    // Member functions

        // Access

            //- Return const access to the element name
            const word& name() const;
            
            //- Return the element location
            vector& position();
            
            //- Return force vector from element
            vector& force();
            
            //- Return relative velocity
            vector& relativeVelocity();
            
            //- Return element angle of attack in degrees
            scalar& angleOfAttack();
            
            //- Return lift coefficient
            scalar& liftCoefficient();
            
            //- Return drag coefficient
            scalar& dragCoefficient();
            
            
        // Manipulation
        
            //- Translate the element's location
            void translate(vector translationVector);
            
            //- Rotate the element about some axis
            void rotate
            (
                vector rotationPoint, 
                vector axis, 
                scalar radians,
                bool rotateVelocity
            );
            void rotate
            (
                vector rotationPoint, 
                vector axis, 
                scalar radians
            );
            
            //- Pitch the element about its half chord
            void pitch(scalar radians);
            
            //- Pitch the element about some fraction of the chord
            void pitch(scalar radians, scalar chordFraction);
            
            //- Set velocity of element
            void setVelocity(vector velocity);
            
            //- Set element speed (velocity magnitude)
            void setSpeed(scalar speed);
            
            //- Set element speed based on rotation
            void setSpeed(vector point, vector axis, scalar omega);
            
            //- Scale the velocity
            void scaleVelocity(scalar scale);
            
            //- Set dynamic stall active
            void setDynamicStallActive(bool active);
            
            //- Set omega for flow curvature correction
            void setOmega(scalar omega);
            

        // Evaluation
            
            //- Calculate forces
            void calculate(vectorField Uin, volVectorField& forceField);
            
            //- Read coefficient data
            void read();
            
            //- Calculate and return moment about specified point
            vector moment(vector point);
            
            
        // Source term addition

            //- Source term to momentum equation
            virtual void addSup
            (
                fvMatrix<vector>& eqn,
                volVectorField& force
            );
            
            //- Add source term to turbulence quantity
            virtual void addTurbulence(fvMatrix<scalar>& eqn, word fieldName);

            //- Source term to compressible momentum equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                volVectorField& force
            );
            
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
