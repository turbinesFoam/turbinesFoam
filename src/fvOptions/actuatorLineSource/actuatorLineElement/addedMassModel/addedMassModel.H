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
    Foam::addedMassModel

Description
    Added mass model based on Strickland (1981) "A vortex model of the
    Darrieus turbine: An analytical and experimental study".

SourceFiles
    addedMassModel.C

\*---------------------------------------------------------------------------*/

#ifndef addedMassModel_H
#define addedMassModel_H

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class addedMassModel Declaration
\*---------------------------------------------------------------------------*/

class addedMassModel
{
    // Private data

        //- Run time reference
        const Time& time_;

        //- Chord length in meters
        scalar chordLength_;

        //- Debug level
        const label& debug;

        //- Previous time value
        scalar timePrev_;

        //- Counter for number of new time steps
        label nNewTimes_;

        //- Angle of attack (radians)
        scalar alpha_;

        //- Previous angle of attack (radians)
        scalar alphaPrev_;

        //- Chordwise component of relative velocity
        scalar chordwiseRelVel_;

        //- Normal component of relative velocity
        scalar normalRelVel_;

        //- Previous value of normal relative velocity
        scalar normalRelVelPrev_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        addedMassModel(const addedMassModel&);

        //- Update previous time step values
        void update();


public:

    // Constructors

        //- Construct from components
        addedMassModel
        (
            const Time& time,
            scalar chordLength,
            const label& debug
        );


    // Selectors

        //- Select from components
        static autoPtr<addedMassModel> New
        (
            const Time& time,
            scalar chordLength,
            const label& debug
        );


    //- Destructor
    ~addedMassModel();


    // Member Functions

        // Access

        // Check

        // Correct

            //- Correct coefficients for added mass effects
            void correct
            (
                scalar& liftCoefficient,
                scalar& dragCoefficient,
                scalar& momentCoefficient,
                scalar alphaRad,
                scalar chordwiseRelVel,
                scalar normalRelVel
            );


        // Write


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
