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

#include "dynamicStallModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(dynamicStallModel, 0);
    defineRunTimeSelectionTable(dynamicStallModel, dictionary);
}
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::fv::dynamicStallModel>
Foam::fv::dynamicStallModel::New
(
    const dictionary& dict,
    const word& modelName,
    const Time& time,
    profileData& profileData
)
{

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn("dynamicStallModel::New(const dictionary&)")
            << "Unknown dynamic stall model type " << modelName
            << nl << nl
            << "Valid model types are :" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<dynamicStallModel>(cstrIter()
    (
        dict,
        modelName,
        time,
        profileData
    ));
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::fv::dynamicStallModel::interpolate
(
    scalar xNew,
    List<scalar>& xOld,
    List<scalar>& yOld
)
{
    if (debug)
    {
        Info<< "dynamicStallModel interpolate function called" << endl;
        Info<< "xNew: " << xNew << endl;
        Info<< "xOld: " << xOld << endl;
        Info<< "yOld: " << yOld << endl;
    }

    label index = 0;
    label indexP = 0;
    label indexM = 0;
    scalar error = 1.0E30;
    forAll(xOld, i)
    {
        scalar diff = mag(xNew - xOld[i]);
        if (diff < error)
        {
            index = i;
            error = diff;
        }
    }
    if (xNew < xOld[index])
    {
        if (index == 0)
        {
            indexP = 1;
            indexM = indexP - 1;
        }
        else
        {
            indexP = index;
            indexM = indexP - 1;
        }
        return yOld[indexM]
               + ((yOld[indexP]
               - yOld[indexM])/(xOld[indexP]
               - xOld[indexM]))*(xNew - xOld[indexM]);
    }
    else if (xNew > xOld[index])
    {
        if (index == xOld.size() - 1)
        {
            indexP = xOld.size() - 1;
            indexM = indexP - 1;
        }
        else
        {
            indexP = index + 1;
            indexM = indexP - 1;
        }
        return yOld[indexM] + ((yOld[indexP]
               - yOld[indexM])/(xOld[indexP]
               - xOld[indexM]))*(xNew - xOld[indexM]);
    }
    else if (xNew == xOld[index])
    {
        return yOld[index];
    }
    else
    {
        return 0.0;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::dynamicStallModel::dynamicStallModel
(
    const dictionary& dict,
    const word& modelName,
    const Time& time,
    profileData& profileData
)
:
    dict_(dict),
    modelName_(modelName),
    time_(time),
    profileData_(profileData),
    coeffs_(dict.subOrEmptyDict(modelName + "Coeffs")),
    startTime_(time.value())
{
    if (debug)
    {
        Info<< modelName << " dynamic stall model created" << endl
            << "    Coeffs:" << endl << coeffs_ << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::dynamicStallModel::~dynamicStallModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::fv::dynamicStallModel::correct(scalar alpha, scalar& cl, scalar& cd)
{}

void Foam::fv::dynamicStallModel::correct
(
    scalar magU,
    scalar alphaDeg,
    scalar& cl,
    scalar& cd,
    scalar& cm
)
{}


void Foam::fv::dynamicStallModel::reduceParallel(bool inMesh){}


// ************************************************************************* //
