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
    Foam::profileData

Description

SourceFiles
    profileDataI.H
    profileData.C

\*---------------------------------------------------------------------------*/

#ifndef profileData_H
#define profileData_H

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;

// Forward declaration of friend functions and operators
class profileData;
Istream& operator>>(Istream&, profileData&);
Ostream& operator<<(Ostream&, const profileData&);


/*---------------------------------------------------------------------------*\
                         Class profileData Declaration
\*---------------------------------------------------------------------------*/

class profileData
{
    // Private data
    
        //- Profile name
        word profileName_;
    
        //- Input dictionary
        dictionary dict_;
    
        //- Reference Reynolds number
        scalar refRe_;
        
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

        //- Disallow default bitwise copy construct
        profileData(const profileData&);


public:

    // Constructors

        //- Construct null
        profileData();

        //- Construct from components
        profileData
        (
            const word& name,
            const dictionary& dict
        );


    // Selectors

        //- Select null constructed
        static autoPtr<profileData> New();


    //- Destructor
    ~profileData();


    // Member Functions

        // Access

        // Check

        // Edit

        // Write


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "profileDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
