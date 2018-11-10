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
    Object to store force and moment coefficient data for a 2-D profile.

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

protected:

    // Protected data

        //- Profile name
        const word name_;

        //- Input dictionary
        const dictionary dict_;

        //- Debug level
        const label& debug;

        //- Table type: "singleRe" (default) or "multiRe"
        const word tableType_;

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

        //- Reynolds number list for multiple Re dataset
        List<scalar> ReList_;

        //- 2-D array of drag coefficient values for multiple Re dataset
        List<List<scalar> > dragCoefficientLists_;

        //- 2-D array of lift coefficient values for multiple Re dataset
        List<List<scalar> > liftCoefficientLists_;

        //- 2-D array of moment coefficient values for multiple Re dataset
        List<List<scalar> > momentCoefficientLists_;

        //- List of static stall angles (deg) for multiple Re dataset
        List<scalar> staticStallAngleList_;

        //- List of zero lift drag coefficients for multiple Re dataset
        List<scalar> zeroLiftDragCoeffList_;

        //- List of zero lift AoAs (deg) for multiple Re dataset
        List<scalar> zeroLiftAngleOfAttackList_;

        //- List of moment coefficients at zero lift for multiple Re dataset
        List<scalar> zeroLiftMomentCoeffList_;

        //- List of normal coefficient slopes (1/rad) for multiple Re dataset
        List<scalar> normalCoeffSlopeList_;

        //- Angle of attack list (deg)
        List<scalar> angleOfAttackList_;

        //- Lift coefficient list at current Re
        List<scalar> liftCoefficientList_;

        //- Drag coefficient list at current Re
        List<scalar> dragCoefficientList_;

        //- Moment coefficient list at current Re
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


    // Protected Member Functions

        //- Read coefficients for single Reynolds number dataset
        void readSingleRe();

        //- Read coefficients for multiple Reynolds number dataset
        void readMultiRe();

        //- Read angle of attack list from first column of a 2-D array
        List<scalar> readAngleOfAttackList(word keyword);

        //- Read 2-D array of coefficient data
        // Data is assumed to have angle of attack values in first column, where
        // the row will be read to xvalues and the column to y values, and data
        // is given as data[y][x]
        void read2DArray
        (
            List<List<scalar> > &data,
            word keyword
        );

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

        //- Set coefficient lists for current Reynolds number by interpolating
        //  the 2-D coefficient arrays
        void interpCoeffLists();

        //- Set properties (static stall angle, etc.) through interpolation
        //  from preconstructed lists
        void interpPropsMultiRe();

        //- Create lists for interpolation with interpPropsMultiRe
        // Only needs to be called once
        void analyzeMultiRe();

        //- Generate sublist from current Re data
        List<scalar> subList
        (
            scalar alphaDegStart,
            scalar alphaDegStop,
            const List<scalar>& fullList
        );

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

            //- Return const reference to dictionary
            const dictionary& dict();

            //- Return const reference to angle of attack list
            const List<scalar>& angleOfAttackList();

            //- Return const reference to lift coefficient list
            const List<scalar>& liftCoefficientList();

            //- Return reference to drag coefficient list
            const List<scalar>& dragCoefficientList();

            //- Return reference to moment coefficient list
            const List<scalar>& momentCoefficientList();

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

            //- Return reference tangential coefficient
            scalar tangentialRefCoefficient
            (
                scalar angleOfAttackDeg,
                scalar inflowVelAngleDeg
            );

            //- Return normal reference coefficient
            scalar normalRefCoefficient
            (
                scalar angleOfAttackDeg,
                scalar inflowVelAngleDeg
            );

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

            //- Convert from lift and drag to tangential reference coefficient
            scalar convertToCRT
            (
                scalar cl,
                scalar cd,
                scalar inflowVelAngleDeg
            );

            //- Convert from lift and drag to normal reference coefficient
            scalar convertToCRN
            (
                scalar cl,
                scalar cd,
                scalar inflowVelAngleDeg
            );


        // Check

            //- Analyze the input data
            void analyze();


        // Edit

            //- Update Reynolds number
            void updateRe(scalar Re);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
