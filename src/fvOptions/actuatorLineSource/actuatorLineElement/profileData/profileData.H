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
        List<scalar> liftCoeffListOrg_;

        //- Unmodified drag coefficient list
        List<scalar> dragCoeffListOrg_;

        //- Unmodified moment coefficient list
        List<scalar> momentCoeffListOrg_;

        //- Angle of attack list (deg)
        List<scalar> angleOfAttackList_;

        //- Lift coefficient list
        List<scalar> liftCoeffList_;

        //- Drag coefficient list
        List<scalar> dragCoeffList_;

        //- Moment coefficient list
        List<scalar> momentCoeffList_;

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

        //- Convert from lift and drag to normal coefficient
        scalar convertToCN
        (
            scalar cl,
            scalar cd,
            scalar angleOfAttackDeg
        );

        //- Convert from lift and drag to chordwise coefficient
        scalar convertToCC
        (
            scalar cl,
            scalar cd,
            scalar angleOfAttackDeg
        );

        //- Convert from normal and chordwise to lift coefficient
        scalar convertToCL
        (
            scalar cn,
            scalar cc,
            scalar angleOfAttackDeg
        );

        //- Convert from normal and chordwise to drag coefficient
        scalar convertToCD
        (
            scalar cn,
            scalar cc,
            scalar angleOfAttackDeg
        );


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

            //- Return const reference to angle of attack list
            const List<scalar>& angleOfAttackList();

            //- Return const reference to lift coefficient list
            const List<scalar>& liftCoeffList();

            //- Return reference to drag coefficient list
            const List<scalar>& dragCoeffList();

            //- Return reference to moment coefficient list
            const List<scalar>& momentCoeffList();

            //- Return subset of angle of attack list
            List<scalar> angleOfAttackList
            (
                scalar alphaDegStart,
                scalar alphaDegStop
            );

            //- Return subset of lift coefficient list
            List<scalar> liftCoeffList
            (
                scalar alphaDegStart,
                scalar alphaDegStop
            );

            //- Return subset of drag coefficient list
            List<scalar> dragCoeffList
            (
                scalar alphaDegStart,
                scalar alphaDegStop
            );

            //- Return subset of moment coefficient list
            List<scalar> momentCoeffList
            (
                scalar alphaDegStart,
                scalar alphaDegStop
            );

            //- Return subset of normal force coefficient list
            List<scalar> normalCoeffList
            (
                scalar alphaDegStart,
                scalar alphaDegStop
            );

            //- Return subset of chordwise force coefficient list
            List<scalar> chordwiseCoeffList
            (
                scalar alphaDegStart,
                scalar alphaDegStop
            );

            //- Lookup lift coefficient
            scalar liftCoeff(scalar angleOfAttackDeg);

            //- Lookup drag coefficient
            scalar dragCoeff(scalar angleOfAttackDeg);

            //- Lookup moment coefficient
            scalar momentCoeff(scalar angleOfAttackDeg);

            //- Return normal force coefficient
            scalar normalCoeff(scalar angleOfAttackDeg);

            //- Return chordwise force coefficient
            scalar chordwiseCoeff(scalar angleOfAttackDeg);

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
