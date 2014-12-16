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

#include "actuatorLineSource.H"
#include "unitConversion.H"
#include "Tuple2.H"
#include "vector.H"
#include "IFstream.H"


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

bool Foam::fv::actuatorLineSource::readFromFile() const
{
    return fName_ != fileName::null;
}


void Foam::fv::actuatorLineSource::interpolateWeights
(
    const scalar& xIn,
    const List<scalar>& values,
    label& i1,
    label& i2,
    scalar& ddx
) const
{
    i2 = 0;
    label nElem = values.size();

    if (nElem == 1)
    {
        i1 = i2;
        ddx = 0.0;
        return;
    }
    else
    {
        while ((values[i2] < xIn) && (i2 < nElem))
        {
            i2++;
        }

        if (i2 == nElem)
        {
            i2 = nElem - 1;
            i1 = i2;
            ddx = 0.0;
            return;
        }
        else
        {
            i1 = i2 - 1;
            ddx = (xIn - values[i1])/(values[i2] - values[i1]);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::actuatorLineSource::actuatorLineSource
(
		const word& name,
		const word& modelType,
		const dictionary& dict,
		const fvMesh& mesh
)
:
    option(name, modelType, dict, mesh),
    profileName_(),
    profileID_(),
    radius_(),
    pitch_(),
    chord_(),
    fName_(fileName::null)
{
    List<Tuple2<word, vector> > data;
    if (readFromFile())
    {
        IFstream is(fName_);
        is  >> data;
    }
    else
    {
        dict.lookup("data") >> data;
    }


    if (data.size() > 0)
    {
        profileName_.setSize(data.size());
        profileID_.setSize(data.size());
        radius_.setSize(data.size());
        pitch_.setSize(data.size());
        chord_.setSize(data.size());

        forAll(data, i)
        {
            profileName_[i] = data[i].first();
            profileID_[i] = -1;
            radius_[i] = data[i].second()[0];
            pitch_[i] = degToRad(data[i].second()[1]);
            chord_[i] = data[i].second()[2];
        }
    }
    else
    {
        FatalErrorIn("Foam::fv::actuatorLineSource::actuatorLineSource(const dictionary&)")
            << "No blade data specified" << exit(FatalError);
    }
}

// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

Foam::fv::actuatorLineSource::~actuatorLineSource()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::List<Foam::word>& Foam::fv::actuatorLineSource::profileName() const
{
    return profileName_;
}


const Foam::List<Foam::label>& Foam::fv::actuatorLineSource::profileID() const
{
    return profileID_;
}


const Foam::List<Foam::scalar>& Foam::fv::actuatorLineSource::radius() const
{
    return radius_;
}


const Foam::List<Foam::scalar>& Foam::fv::actuatorLineSource::pitch() const
{
    return pitch_;
}


const Foam::List<Foam::scalar>& Foam::fv::actuatorLineSource::chord() const
{
    return chord_;
}


Foam::List<Foam::label>& Foam::fv::actuatorLineSource::profileID()
{
    return profileID_;
}


void Foam::fv::actuatorLineSource::interpolate
(
    const scalar radius,
    scalar& pitch,
    scalar& chord,
    label& i1,
    label& i2,
    scalar& invDr
) const
{
    interpolateWeights(radius, radius_, i1, i2, invDr);

    pitch = invDr*(pitch_[i2] - pitch_[i1]) + pitch_[i1];
    chord = invDr*(chord_[i2] - chord_[i1]) + chord_[i1];
}


// ************************************************************************* //
