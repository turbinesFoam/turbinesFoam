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
    Foam::fv::actuatorLineSource

Description
    Actuator line class, which is a collection of actuator line elements.

SourceFiles
    actuatorLineSource.C

\*---------------------------------------------------------------------------*/

#ifndef actuatorLineSource_H
#define actuatorLineSource_H

#include "List.H"
#include "dictionary.H"
#include "vector.H"
#include "actuatorLineElement.H"
#include "fvOption.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                         Class actuatorLineSource Declaration
\*---------------------------------------------------------------------------*/

class actuatorLineSource
:
    public option
{

protected:

    // Protected data
        
        //- Total length of geometry in meters
        scalar totalLength_;
        
        //- Tip effect
        scalar tipEffect_;
        
        //- Sectional coefficient data for the foil profile
        List<List<scalar> > coefficientData_;
        
        //- Element geometry
        // point, spanDirection, chordLength, chordRefDir, pitch 
        List<List<List<scalar> > > elementGeometry_;
        
        //- Number of blade elements
        label nElements_;
        
        //- Free stream velocity
        vector freeStreamVelocity_;
        
        //- Free stream direction
        vector freeStreamDirection_;
        
        //- Total force vector from all elements
        vector force_;
        
        //- Force field from all elements
        volVectorField forceField_;
        
        //- List of actuator line elements
        PtrList<actuatorLineElement> elements_;


    // Protected Member Functions
        
        //- Create actuator line elements
        void createElements();
        
        //- Read dictionary
        bool read(const dictionary& dict);


public:

    //- Runtime type information
    TypeName("actuatorLineSource");

    // Selectors

    //- Return a reference to the selected fvOption model
    static autoPtr<actuatorLineSource> New
    (
        const word& name,
        const dictionary& dict,
        const fvMesh& mesh
    );
        

    //- Constructor
    actuatorLineSource
	(
		const word& name,
		const word& modelType,
		const dictionary& dict,
		const fvMesh& mesh
	);


    //- Destructor
    virtual ~actuatorLineSource();


    // Member functions

        // Access
            
            //- Return the total force vector
            vector& force();
            
            //- Return the force field
            volVectorField& forceField();


        // Edit
            
            //- Translate the actuator line
            void translate(vector translation);
            
            //- Rotate the actuator line about a specified axis
            void rotate(vector rotationPoint, vector axis, scalar radians);
            
            //- Pitch the blade about its half chord
            void pitch(scalar radians);


        // Evaluation
        
            //- Compute the moment about a given point
            vector moment(vector point);
            
        // IO
        
            //- Print dictionary values
            virtual void printCoeffs() const;


        // Source term addition

            //- Source term to momentum equation
            virtual void addSup
            (
                fvMatrix<vector>& eqn,
                const label fieldI
            );
            
            //- Source term to turbulence scalars
            virtual void addSup
            (
                fvMatrix<scalar>& eqn,
                const label fieldI
            );

            //- Source term to compressible momentum equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                const label fieldI
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
