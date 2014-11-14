/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::actuatorLineElementList

Description
    Base class for profile models

SourceFiles
    actuatorLineElementList.C

\*---------------------------------------------------------------------------*/

#ifndef actuatorLineElementList_H
#define actuatorLineElementList_H

#include "PtrList.H"
#include "actuatorLineElement.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class actuatorLineElementList Declaration
\*---------------------------------------------------------------------------*/

class actuatorLineElementList
:
    public PtrList<actuatorLineElement>
{

protected:

    // Protected data

        //- Dictionary
        const dictionary dict_;


public:

    //- Constructor
    actuatorLineElementList(const dictionary& dict, const bool readFields = true);

    //- Destructor
    ~actuatorLineElementList();


    // Member Functions

        //- Set blade->profile addressing
        void connectBlades
        (
            const List<word>& names,
            List<label>& addr
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
