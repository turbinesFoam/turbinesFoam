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

\*---------------------------------------------------------------------------*/

#include "interpolateUtils.H"

label Foam::interpolateUtils::binarySearch
(
     const List<scalar>& list,
     const scalar value
)
{
    // Find the closest index with a list value below value using a binary
    // search algorithm, which should find the value in O(ln(N)) time
    // suitable for large lists
    label index = 0;
    label listsize = list.size();
    for
    (
        label currentsize = listsize >> 1;
        currentsize > 1 && index + currentsize < listsize;
        currentsize = ((currentsize+1)>>1)
    )
    {
        if (value > list[index + currentsize])
        {
            index += currentsize;
        }
    }

    // The last case has to be run separately, as currentsize = 1 cannot be
    // handled in the loop due to the currentsize decrement
    if (index + 1 < listsize && value > list[index+1])
    {
        index += 1;
    }
    return index;
}

label Foam::interpolateUtils::linearSearch
(
    const List<scalar>& list,
    const scalar value,
    label startvalue
)
{
    // Find the closest index with a list value below the value using a linear
    // search algorithm, which should find the value in O(N) time
    // Suitable for small lists, or if a good startvalue is available
    // result should be indentical to binarySearch
    label listsize = list.size();
    if (startvalue >= listsize)
    {
        startvalue = listsize - 1;
    }
    if (list[startvalue] < value)
    {
        while (startvalue + 1 < listsize && list[startvalue+1] < value)
        {
            startvalue++;
        }
    }
    else
    {
        while (startvalue > 0 && list[startvalue] >= value)
        {
            startvalue--;
        }
    }
    return startvalue;
}

scalar Foam::interpolateUtils::getPart
(
    const scalar xNew,
    const List<scalar>& xList,
    label &xIndex
)
{
    scalar xPart;

    // If the value is outside the interpolation region, the edge value will be
    // used

    // If value is smaller than lowest value, use the lowest value
    if (xIndex == 0 && xNew < xList[0])
    {
        xPart = 0;
    }
    // If it is larger than largest value
    else if (xIndex + 1 == xList.size())
    {
        // Decrease value one step, but only use the final value, i.e., the
        // highest value
        xIndex--;
        xPart = 1;
    }
    else
    {
        xPart = (xNew - xList[xIndex])/(xList[xIndex+1] - xList[xIndex]);
    }

    return xPart;
}


scalar Foam::interpolateUtils::interpolate1D
(
    const scalar xNew,
    const List<scalar>& xList,
    const List<scalar>& data,
    label xIndex
)
{
    // Index is known
    return interpolate1D
    (
        getPart(xNew, xList, xIndex),
        data,
        xIndex
    );
}


scalar Foam::interpolateUtils::interpolate1D
(
    const scalar xPart,
    const List<scalar>& data,
    label xIndex
)
{
    // Interpolate fraction is known
    return (data[xIndex]*(1 - xPart) + data[xIndex + 1]*xPart);
}


scalar Foam::interpolateUtils::interpolate1D
(
    const scalar xNew,
    const List<scalar>& xList,
    const List<scalar>& data
)
{
    // General 1-D interpolation
    return interpolate1D
    (
        xNew,
        xList,
        data,
        binarySearch(xList, xNew)
    );
}


scalar Foam::interpolateUtils::interpolate2D
(
    const scalar xNew,
    const scalar yNew,
    const List<scalar>& xList,
    const List<scalar>& yList,
    const List<List<scalar> >& data,
    label xIndex,
    label yIndex
)
{
    // Index values are known
    return interpolate2D
    (
        getPart(xNew, xList, xIndex),
        getPart(yNew, yList, yIndex),
        data,
        xIndex,
        yIndex
    );
}


scalar Foam::interpolateUtils::interpolate2D
(
    const scalar xPart,
    const scalar yPart,
    const List<List<scalar> >& data,
    label xIndex,
    label yIndex
)
{
    // Interpolation fractions have been calculated
    return (data[yIndex][xIndex]*(1 - xPart) +
            data[yIndex][xIndex+1]*xPart)*(1 - yPart) +
           (data[yIndex+1][xIndex]*(1 - xPart) +
            data[yIndex + 1][xIndex + 1]*xPart)*yPart;
}


scalar Foam::interpolateUtils::interpolate2D
(
    const scalar xNew,
    const scalar yNew,
    const List<scalar>& xList,
    const List<scalar>& yList,
    const List<List<scalar> >& data
)
{
    return interpolate2D
    (
        xNew,
        yNew,
        xList,
        yList,
        data,
        binarySearch(xList, xNew),
        binarySearch(yList, yNew)
    );
}


// ************************************************************************* //
