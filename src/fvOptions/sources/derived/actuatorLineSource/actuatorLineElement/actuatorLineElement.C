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

\*---------------------------------------------------------------------------*/

#include "actuatorLineElement.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMatrices.H"
#include "geometricOneField.H"
#include "syncTools.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(actuatorLineElement, 0);
    defineRunTimeSelectionTable(actuatorLineElement, dictionary);
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

bool Foam::fv::actuatorLineElement::readFromFile() const
{
    return fName_ != fileName::null;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::actuatorLineElement::actuatorLineElement
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    option(name, modelType, dict, mesh),
    dict_(dict),
    name_(name),
    fName_(fileName::null)
{
    dict.readIfPresent("fileName", fName_);
}

// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

Foam::fv::actuatorLineElement::~actuatorLineElement()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::word& Foam::fv::actuatorLineElement::name() const
{
    return name_;
}


void Foam::fv::actuatorLineElement::calculate()
{
	Info<< "Calculating force contribution from actuator line element" << endl;
	// Calculate local wind velocity
	// Calculate relative velocity
	// Calculate angle of attack
	// Lookup lift and drag coefficients
	// Calculate force
	// Calculate force field 
}


Foam::vector& Foam::fv::actuatorLineElement::force()
{
    calculate();
    return force_;
}


void Foam::fv::actuatorLineElement::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    // Read the reference density for incompressible flow
    coeffs_.lookup("rhoRef") >> rhoRef_;

    calculate();

    // Add source to rhs of eqn
    //eqn -= forceField_;
}


void Foam::fv::actuatorLineElement::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    calculate();

    // Add source to rhs of eqn
    //eqn -= forceField_;
}


Foam::autoPtr<Foam::fv::actuatorLineElement> Foam::fv::actuatorLineElement::New
(
    const dictionary& dict
)
{
    const word& modelName(dict.dictName());

    const word modelType(dict.lookup("type"));

    Info<< "    - creating " << modelType << " profile " << modelName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn("actuatorLineElement::New(const dictionary&)")
            << "Unknown profile model type " << modelType
            << nl << nl
            << "Valid model types are :" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<actuatorLineElement>(cstrIter()(dict, modelName));
}


// ************************************************************************* //
