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

        //- Corresponding profile name per section
        List<word> profileName_;

        //- Corresponding profile ID per section
        List<label> profileID_;

        //- List of blade element axial distance from origin
        List<scalar> axialDistance_;
        
        //- List of blade element radii in meters
        List<scalar> radius_;

        //- List of blade element initial azimuthal positions in degrees
        List<scalar> initialAzimuth_;

        //- Chord mounting location relative to radius, as a fraction of chord
        List<scalar> chordMount_;

        //- Preset blade pitch in degrees, eventually converted to radians
        List<scalar> pitch_;

        //- Chord in meters
        List<scalar> chord_;

        //- File name (optional)
        fileName fName_;
        
        //- List of actuator line elements
        PtrList<actuatorLineElement> actuatorLineElements_;
        
        //- Tip effect
        scalar tipEffect_;
        
        //- Sectional coefficient data for the foil profile
        List<List<scalar> > coefficientData_;
        
        //- Element geometry
        // point, spanDirection, chordLength, pitch (from freeStreamVelocity_)
        List<List<List<scalar> > > elementGeometry_;
        
        //- Number of blade elements
        label nElements_;
        
        //- Free stream velocity
        vector freeStreamVelocity_;
        
        //- List of actuator line elements
        PtrList<actuatorLineElement> elements_;


    // Protected Member Functions

        //- Return the interpolation indices and gradient
        void interpolateWeights
        (
            const scalar& xIn,
            const List<scalar>& values,
            label& i1,
            label& i2,
            scalar& ddx
        ) const;
        
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

            //- Return const access to the profile name list
            const List<word>& profileName() const;

            //- Return const access to the profile ID list
            const List<label>& profileID() const;

            //- Return const access to the radius list
            const List<scalar>& radius() const;

            //- Return const access to the pitch list
            const List<scalar>& pitch() const;
            
            //- Return const access to the initialAzimuth list
            const List<scalar>& initialAzimuth() const;

            //- Return const access to the chord list
            const List<scalar>& chord() const;
            
            //- Return const access to the chordMount list
            const List<scalar>& chordMount() const;


        // Edit

            //- Return non-const access to the profile ID list
            List<label>& profileID();
            
            //- Translate the actuator line
            void translate(vector translation);
            
            //- Rotate the actuator line about a specified axis
            void rotate(vector rotationPoint, vector axis, scalar radians);


        // Evaluation

            //- Return the twist and chord for a given radius
            virtual void interpolate
            (
                const scalar radius,
                scalar& twist,
                scalar& chord,
                label& i1,
                label& i2,
                scalar& invDr
            ) const;
            
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
