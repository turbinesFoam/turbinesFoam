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
    Foam::profileData

Description

SourceFiles
    profileData.C

\*---------------------------------------------------------------------------*/

#ifndef profileData_H
#define profileData_H

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class profileData Declaration
\*---------------------------------------------------------------------------*/

class profileData
{
    // Private data

        //- Profile name
        const word name_;

        //- Input dictionary
        const dictionary dict_;

        //- Debug level
        const label& debug;

        //- Specified Reynolds number
        scalar Re_;

        //- Reference Reynolds number
        scalar ReRef_;

        //- Unmodified angle of attack list (deg)
        List<scalar> angleOfAttackListOrg_;

        //- Unmodified lift coefficient list
        List<scalar> liftCoefficientListOrg_;

        //- Unmodified drag coefficient list
        List<scalar> dragCoefficientListOrg_;

        //- Unmodified moment coefficient list
        List<scalar> momentCoefficientListOrg_;

        //- Angle of attack list (deg)
        List<scalar> angleOfAttackList_;

        //- Lift coefficient list
        List<scalar> liftCoefficientList_;

        //- Drag coefficient list
        List<scalar> dragCoefficientList_;

        //- Moment coefficient list
        List<scalar> momentCoefficientList_;

        //- Switch for Reynolds number corrections
        bool correctRe_;

        //- Static stall angle (deg)
        scalar staticStallAngle_;

        //- Drag coefficient at zero lift
        scalar zeroLiftDragCoeff_;

        //- Angle of attack at zero lift (deg)
        scalar zeroLiftAngleOfAttack_;

        //- Moment coefficient at zero lift
        scalar zeroLiftMomentCoeff_;

        //- Slope of normal force coefficient (1/rad)
        scalar normalCoeffSlope_;


    // Private Member Functions

        //- Read dictionary
        void read();

        //- Interpolate a scalar value
        scalar interpolate
        (
            scalar xNew,
            List<scalar>& xOld,
            List<scalar>& yOld
        );

        //- Calculate the static stall angle
        void calcStaticStallAngle();

        //- Calculate the zero lift drag coefficient
        void calcZeroLiftDragCoeff();

        //- Calculate zero lift angle of attack
        void calcZeroLiftAngleOfAttack();

        //- Calculate zero lift moment coefficient
        void calcZeroLiftMomentCoeff();

        //- Calculate normal force coefficient slope
        void calcNormalCoeffSlope();

        //- Disallow default bitwise copy construct
        profileData(const profileData&);


public:

    // Constructors

        //- Construct from components
        profileData
        (
            const word& name,
            const dictionary& dict,
            const label& debug
        );


    // Selectors

        //- Select from components
        static autoPtr<profileData> New
        (
            const word& name,
            const dictionary& dict,
            const label& debug
        );


    //- Destructor
    ~profileData();


    // Member Functions

        // Access

            //- Return const reference to dictionary
            const dictionary& dict();

            //- Return reference to angle of attack list
            List<scalar>& angleOfAttackList();

            //- Return reference to lift coefficient list
            List<scalar>& liftCoefficientList();

            //- Return reference to drag coefficient list
            List<scalar>& dragCoefficientList();

            //- Return reference to moment coefficient list
            List<scalar>& momentCoefficientList();

            //- Return subset of angle of attack list
            List<scalar> angleOfAttackList
            (
                scalar alphaDegStart,
                scalar alphaDegStop
            );

            //- Return subset of lift coefficient list
            List<scalar> liftCoefficientList
            (
                scalar alphaDegStart,
                scalar alphaDegStop
            );

            //- Return subset of drag coefficient list
            List<scalar> dragCoefficientList
            (
                scalar alphaDegStart,
                scalar alphaDegStop
            );

            //- Return subset of moment coefficient list
            List<scalar> momentCoefficientList
            (
                scalar alphaDegStart,
                scalar alphaDegStop
            );

            //- Return subset of normal force coefficient list
            List<scalar> normalCoefficientList
            (
                scalar alphaDegStart,
                scalar alphaDegStop
            );

            //- Return subset of chordwise force coefficient list
            List<scalar> chordwiseCoefficientList
            (
                scalar alphaDegStart,
                scalar alphaDegStop
            );

            //- Lookup lift coefficient
            scalar liftCoefficient(scalar angleOfAttackDeg);

            //- Lookup drag coefficient
            scalar dragCoefficient(scalar angleOfAttackDeg);

            //- Lookup moment coefficient
            scalar momentCoefficient(scalar angleOfAttackDeg);

            //- Return normal force coefficient
            scalar normalCoefficient(scalar angleOfAttackDeg);

            //- Return chordwise force coefficient
            scalar chordwiseCoefficient(scalar angleOfAttackDeg);

            //- Return static stall angle in radians
            scalar staticStallAngleRad();

            //- Return zero lift drag coefficient
            scalar zeroLiftDragCoeff();

            //- Return zero lift angle of attack (deg)
            scalar zeroLiftAngleOfAttack();

            //- Return zero lift moment coefficient
            scalar zeroLiftMomentCoeff();

            //- Return normal coefficient slope (1/rad)
            scalar normalCoeffSlope();

            //- Return Reynolds number
            scalar Re();

            //- Indicate if Reynolds number correction is active
            bool correctRe();


        // Check

            //- Analyze the input data
            void analyze();


        // Edit

            //- Update Reynolds number
            void updateRe(scalar Re);


        // Write


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
