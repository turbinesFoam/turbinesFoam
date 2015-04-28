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
    Foam::fv::dynamicStallModel

Description

SourceFiles
    dynamicStallModel.C

\*---------------------------------------------------------------------------*/

#ifndef dynamicStallModel_H
#define dynamicStallModel_H

#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "dictionary.H"
#include "vector.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                         Class dynamicStallModel Declaration
\*---------------------------------------------------------------------------*/

class dynamicStallModel
{
    
protected:

    // Protected data
    
        //- Model name
        word modelName_;


public:

    //- Runtime type information
    TypeName("dynamicStallModel");

        // Declare run-time new constructor selection table
        declareRunTimeNewSelectionTable
        (
            autoPtr,
            dynamicStallModel,
            dictionary,
            (
                const dictionary& dict,
                const word& modelName
            ),
            (dict, modelName)
        );

    // Constructors

        //- Construct from components
        dynamicStallModel(const dictionary&, const word& modelName);

        //- Construct as copy
        dynamicStallModel(const dynamicStallModel&);


    // Selectors

        //- Select null constructed
        static autoPtr<dynamicStallModel> New
        (
            const dictionary&, 
            const word& modelName
        );


    //- Destructor
    virtual ~dynamicStallModel();


    // Member Functions

        // Access

        // Check

        // Edit
        
        // Evaluation
        
            // Correct lift and drag coefficients
            virtual void correct(scalar alphaDeg, scalar& cl, scalar& cd);

        // Write


    // Member Operators

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
