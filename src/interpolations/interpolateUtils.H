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

Namespace
    Foam::interpolateUtils

Description
    Useful functions for interpolating data from input files.

SourceFiles
    interpolateUtils.C

\*---------------------------------------------------------------------------*/

#ifndef interpolateUtils_H
#define interpolateUtils_H

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace interpolateUtils
{

    //- Find index for interpolation
    label binarySearch
    (
        const List<scalar>& list,
        const scalar value
    );

    //- Find index for interpolation
    label linearSearch
    (
        const List<scalar>& list,
        const scalar value,
        label startvalue = 0
    );

    //- Get interpolation fraction. Note: Does not work if any of the lists
    //  only contain one variable.
    scalar getPart
    (
        const scalar xNew,
        const List<scalar>& xList,
        label &xIndex
    );


    //- 1-D iterpolation functions

        //- Interpolate from 1-D data with known index
        scalar interpolate1D
        (
            const scalar xNew,
            const List<scalar>& xList,
            const List<scalar>& data,
            label xIndex
        );

        //- Interpolate from 1-D data with known fraction
        scalar interpolate1D
        (
            const scalar xPart,
            const List<scalar>& data,
            label xIndex
        );

        //- Interpolate from 1-D data
        scalar interpolate1D
        (
            const scalar xNew,
            const List<scalar>& xList,
            const List<scalar>& data
        );


    //- 2-D Interpolation functions
    //  Note: Will not work if any of the lists only contain one variable.
    //  The format of data should be data[y][x]

        //- Interpolate from 2-D data with known indices
        scalar interpolate2D
        (
            const scalar xNew,
            const scalar yNew,
            const List<scalar>& xList,
            const List<scalar>& yList,
            const List<List<scalar> >& data,
            label xIndex,
            label yIndex
        );

        //- Interpolate from 2-D data with known fractions
        scalar interpolate2D
        (
            const scalar xPart,
            const scalar yPart,
            const List<List<scalar> >& data,
            label xIndex,
            label yIndex
        );

        //- Interpolate from 2-D data
        scalar interpolate2D
        (
            const scalar xNew,
            const scalar yNew,
            const List<scalar>& xList,
            const List<scalar>& yList,
            const List<List<scalar> >& data
        );
} // End namespace interpolateUtils


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
