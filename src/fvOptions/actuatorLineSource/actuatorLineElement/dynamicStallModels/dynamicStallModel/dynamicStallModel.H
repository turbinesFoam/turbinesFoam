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
    Foam::fv::dynamicStallModel

Description
    Generic base class for dynamic stall model.

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
#include "profileData.H"

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

        //- Dictionary
        const dictionary dict_;

        //- Model name
        const word modelName_;

        //- Runtime reference
        const Time& time_;

        //- Profile data reference
        profileData& profileData_;

        //- Coefficients subdict
        const dictionary coeffs_;

        //- Starting time for the model
        scalar startTime_;


    // Protected member functions

        //- Interpolate a scalar value
        scalar interpolate
        (
            scalar xNew,
            List<scalar>& xOld,
            List<scalar>& yOld
        );


public:

    //- Runtime type information
    TypeName("dynamicStallModel");

        //- Declare run-time constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            dynamicStallModel,
            dictionary,
            (
                const dictionary& dict,
                const word& modelName,
                const Time& time,
                profileData& profileData
            ),
            (dict, modelName, time, profileData)
        );

    // Constructors

        //- Construct from components
        dynamicStallModel
        (
            const dictionary&,
            const word& modelName,
            const Time& time,
            profileData& profileData
        );


    // Selectors

        //- Select from components
        static autoPtr<dynamicStallModel> New
        (
            const dictionary& dict,
            const word& modelName,
            const Time& time,
            profileData& profileData
        );


    //- Destructor
    virtual ~dynamicStallModel();


    // Member Functions

        // Access

        // Check

        // Edit

        // Evaluation

            //- Correct lift and drag coefficients
            virtual void correct(scalar alphaDeg, scalar& cl, scalar& cd);

            //- Correct lift and drag coefficients
            virtual void correct
            (
                scalar magU,
                scalar alphaDeg,
                scalar& cl,
                scalar& cd,
                scalar& cm
            );

        // Write

        // Parallel running

            //- Reduce to set values equal on all processors
            virtual void reduceParallel(bool inMesh);


    // Member Operators

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
