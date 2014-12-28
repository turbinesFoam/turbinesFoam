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
#include "addToRunTimeSelectionTable.H"
#include "vector.H"
#include "fvMatrices.H"
#include "geometricOneField.H"
#include "syncTools.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(actuatorLineSource, 0);
    addToRunTimeSelectionTable
    (
        option,
        actuatorLineSource,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //


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


Foam::scalar Foam::fv::actuatorLineSource::interpolate
(
    scalar xNew, 
    List<scalar>& xOld, 
    List<scalar>& yOld
)
{
    label index = 0;
    label indexP = 0;
    label indexM = 0;
    scalar error = 1.0E30;
    forAll(xOld, i)
    {
        scalar diff = mag(xNew - xOld[i]);
        if(diff < error)
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
    read(dict_);
    createElements();
}

// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

Foam::fv::actuatorLineSource::~actuatorLineSource()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::List<Foam::word>& Foam::fv::actuatorLineSource::profileName() const
{
    return profileName_;
}


void Foam::fv::actuatorLineSource::printCoeffs() const
{
    // Print turbine properties
    Info<< "Actuator line properties:" << endl;
    Info<< "Coefficient data:" << endl;
    Info<< coefficientData_ << endl;
    Info<< "First item of element geometry:" << endl;
    Info<< elementGeometry_[0] << endl;
}


bool Foam::fv::actuatorLineSource::read(const dictionary& dict)
{
    if (option::read(dict))
    {
        
        coeffs_.lookup("fieldNames") >> fieldNames_;
        applied_.setSize(fieldNames_.size(), false);

        // Get foil information
        coeffs_.lookup("coefficientData") >> coefficientData_;
        coeffs_.lookup("tipEffect") >> tipEffect_;
        coeffs_.lookup("elementGeometry") >> elementGeometry_;
        coeffs_.lookup("nElements") >> nElements_;
        
        coeffs_.lookup("freeStreamVelocity") >> freeStreamVelocity_;
        
        if (debug)
        {
            Info<< "Debugging for actuatorLineSource on" << endl;
            printCoeffs();
        }

        return true;
    }
    else
    {
        return false;
    }
}


void Foam::fv::actuatorLineSource::createElements()
{
	elements_.setSize(nElements_);
    
    label nGeometryPoints = elementGeometry_.size();
    List<vector> points(nGeometryPoints);
    List<vector> spanDirs(nGeometryPoints);
    List<scalar> chordLengths(nGeometryPoints);
    List<scalar> pitches(nGeometryPoints);
    
    for (int i = 0; i < nGeometryPoints; i++)
    {
        // Extract geometry point
        scalar x = elementGeometry_[i][0][0];
        scalar y = elementGeometry_[i][0][1];
        scalar z = elementGeometry_[i][0][2];
        points[i] = vector(x, y, z);
        // Read span direction
        x = elementGeometry_[i][1][0];
        y = elementGeometry_[i][1][1];
        z = elementGeometry_[i][1][2];
        spanDirs[i] = vector(x, y, z);
        // Read chord length
        chordLengths[i] = elementGeometry_[i][2][0];
        // Read pitch
        pitches[i] = elementGeometry_[i][3][0];
    }
    
    if (debug)
    {
        Info<< "Points:" << endl << points << endl;
        Info<< "Span directions:" << endl << spanDirs << endl;
        Info<< "Chord lengths:" << endl << chordLengths << endl;
        Info<< "Pitches:" << endl << pitches << endl;
    }
	
    for (label i = 0; i < nElements_; i++)
    {
        std::stringstream ss;
        ss << i;
        string str = ss.str();
        const word name = name_ + "Element" + str;

        // Sample values -- should be calculated from elementGeometry
        vector position;
        scalar chordLength;
        vector chordDirection;
        scalar spanLength;
        vector spanDirection;
        scalar pitch;
        
        if (i == 0)
        {
            position = points[i];
            chordLength = chordLengths[i];
            spanLength = mag(points[i] - points[i+1]);
            spanDirection = spanDirs[i];
            pitch = pitches[i];
        }
        else if (i == nElements_ - 1)
        {
            label ind = nGeometryPoints - 1;
            position = points[ind];
            chordLength = chordLengths[ind];
            spanLength = mag(points[i] - points[i-1]);
            spanDirection = spanDirs[ind];
            pitch = pitches[ind];
        }
        else
        {
            position = vector(0, 0, 0);
            chordLength = 0.0;
            spanLength = 0.0;
            spanDirection = vector(0, 0, 0);
            pitch = 0.0;
        }
        
        // Chord direction points into free stream then rotated by pitch
        chordDirection = -freeStreamVelocity_;
        
        // Create a dictionary for this actuatorLineElement
        dictionary dict;
        dict.add("coefficientData", coefficientData_);
        dict.add("chordLength", chordLength);
        dict.add("chordDirection", chordDirection);
        dict.add("spanLength", spanLength);
        dict.add("spanDirection", spanDirection);
        
        if (debug)
        {
            Info<< "Creating actuatorLineElement: " << name << endl;
            Info<< "Position: " << position << endl;
            Info<< "Chord length: " << chordLength << endl;
            Info<< "Chord direction: " << chordDirection << endl;
            Info<< "Span length: " << spanLength << endl;
            Info<< "Span direction: " << spanDirection << endl;
        }
        
        actuatorLineElement* element = new actuatorLineElement(name, dict, mesh_);
        elements_.set(i, element);
        pitch = pitch/180.0*Foam::constant::mathematical::pi;
        elements_[i].pitch(pitch);
    }
}


void Foam::fv::actuatorLineSource::rotate(vector rotationPoint, vector axis,
                                          scalar radians)
{
    forAll(elements_, i)
    {
        elements_[i].rotate(rotationPoint, axis, radians);
    }
} 


void Foam::fv::actuatorLineSource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    volVectorField force
    (
        IOobject
        (
            name_ + ":force",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedVector
        (
            "zero",
            eqn.dimensions()/dimVolume,
            vector::zero
        )
    );

    // Read the reference density for incompressible flow
    //coeffs_.lookup("rhoRef") >> rhoRef_;

    //const vectorField Uin(inflowVelocity(eqn.psi()));
    //calculate();
    
    for (int i = 0; i < nElements_; i++)
    {
        elements_[i].addSup(eqn, force);
    }

    // Add source to rhs of eqn
    eqn -= force;

    if (mesh_.time().outputTime())
    {
        force.write();
    }
    
    Info<< "Force contribution from actuator line " << name_ << ":" << endl;
    Info<< force.average() << endl;
}


void Foam::fv::actuatorLineSource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    volVectorField force
    (
        IOobject
        (
            name_ + ":force",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedVector
        (
            "zero",
            eqn.dimensions()/dimVolume,
            vector::zero
        )
    );

    //const vectorField Uin(inflowVelocity(eqn.psi()));
    //calculate();

    // Add source to rhs of eqn
    eqn -= force;

    if (mesh_.time().outputTime())
    {
        force.write();
    }
}

// ************************************************************************* //
