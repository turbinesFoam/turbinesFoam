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
    forceField_
    (
        IOobject
        (
            "force." + name_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector
        (
            "force", 
            dimForce/dimVolume/dimDensity, 
            vector::zero
        )
    )
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

        // Look up information in dictionary
        coeffs_.lookup("coefficientData") >> coefficientData_;
        coeffs_.lookup("tipEffect") >> tipEffect_;
        coeffs_.lookup("elementGeometry") >> elementGeometry_;
        coeffs_.lookup("nElements") >> nElements_;
        coeffs_.lookup("freeStreamVelocity") >> freeStreamVelocity_;
        freeStreamDirection_ = freeStreamVelocity_/mag(freeStreamVelocity_);
        
        
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
    label nGeometrySegments = nGeometryPoints - 1;
    label nElementsPerSegment = nElements_/nGeometrySegments;
    if (nElements_ % nGeometrySegments)
    {
        // Need to have integer number of elements per geometry segment
        FatalErrorIn("void actuatorLineSource::createElements()")
            << "Number of actuator line elements must be multiple of the number of actuator line geometry segments" 
            << abort(FatalError);
    }
    List<vector> points(nGeometryPoints);
    List<vector> spanDirs(nGeometryPoints);
    List<scalar> chordLengths(nGeometryPoints);
    List<vector> chordRefDirs(nGeometryPoints);
    List<scalar> pitches(nGeometryPoints);
    totalLength_ = 0.0;
    
    for (int i = 0; i < nGeometryPoints; i++)
    {
        // Extract geometry point
        scalar x = elementGeometry_[i][0][0];
        scalar y = elementGeometry_[i][0][1];
        scalar z = elementGeometry_[i][0][2];
        points[i] = vector(x, y, z);
        if (i > 0) totalLength_ += mag(points[i] - points[i-1]);
        // Read span direction
        x = elementGeometry_[i][1][0];
        y = elementGeometry_[i][1][1];
        z = elementGeometry_[i][1][2];
        spanDirs[i] = vector(x, y, z);
        // Read chord length
        chordLengths[i] = elementGeometry_[i][2][0];
        // Read chord ref dir
        x = elementGeometry_[i][3][0];
        y = elementGeometry_[i][3][1];
        z = elementGeometry_[i][3][2];
        chordRefDirs[i] = vector(x, y, z);
        // Read pitch
        pitches[i] = elementGeometry_[i][4][0];
    }
    
    // Calculate span length of each element
    scalar spanLength = totalLength_/nElements_;
    
    // Lookup initial element velocities if present
    List<vector> initialVelocities(nElements_, vector::zero);
    coeffs_.readIfPresent("initialVelocities", initialVelocities);
    
    if (debug)
    {
        Info<< "Total length: " << totalLength_ << endl;
        Info<< "Element span length: " << spanLength << endl;
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

        // Actuator point geometry to be calculated from elementGeometry
        label geometrySegmentIndex = i/nElementsPerSegment;
        label pointIndex = i % nElementsPerSegment;
        vector position;
        scalar chordLength;
        vector chordDirection;
        vector spanDirection;
        scalar pitch;
        vector initialVelocity;
        
        // Linearly interpolate position
        vector point1 = points[geometrySegmentIndex];
        vector point2 = points[geometrySegmentIndex + 1];
        vector segment = point2 - point1;
        position = point1 
                 + segment/nElementsPerSegment*pointIndex
                 + segment/nElementsPerSegment/2;
                 
        // Linearly interpolate chordLength
        scalar chordLength1 = chordLengths[geometrySegmentIndex];
        scalar chordLength2 = chordLengths[geometrySegmentIndex + 1];
        scalar deltaChordTotal = chordLength2 - chordLength1;
        chordLength = chordLength1 
                    + deltaChordTotal/nElementsPerSegment*pointIndex
                    + deltaChordTotal/nElementsPerSegment/2;
                    
        // Linearly interpolate spanDirection
        vector spanDir1 = spanDirs[geometrySegmentIndex];
        vector spanDir2 = spanDirs[geometrySegmentIndex + 1];
        vector deltaSpanTotal = spanDir2 - spanDir1;
        spanDirection = spanDir1 
                      + deltaSpanTotal/nElementsPerSegment*pointIndex
                      + deltaSpanTotal/nElementsPerSegment/2;

        // Linearly interpolate section pitch
        scalar pitch1 = pitches[geometrySegmentIndex];
        scalar pitch2 = pitches[geometrySegmentIndex + 1];
        scalar deltaPitchTotal = pitch2 - pitch1;
        pitch = pitch1 
              + deltaPitchTotal/nElementsPerSegment*pointIndex
              + deltaPitchTotal/nElementsPerSegment/2;
              
        // Linearly interpolate element velocity
        vector vel1 = initialVelocities[geometrySegmentIndex];
        vector vel2 = initialVelocities[geometrySegmentIndex + 1];
        vector deltaVelTotal = vel2 - vel1;
        initialVelocity = vel1 
                        + deltaVelTotal/nElementsPerSegment*pointIndex
                        + deltaVelTotal/nElementsPerSegment/2;
        
        // Linearly interpolate chordDirection
        vector chordDir1 = chordRefDirs[geometrySegmentIndex];
        vector chordDir2 = chordRefDirs[geometrySegmentIndex + 1];
        vector deltaChordDirTotal = chordDir2 - chordDir1;
        chordDirection = chordDir1 
                       + deltaChordDirTotal/nElementsPerSegment*pointIndex
                       + deltaChordDirTotal/nElementsPerSegment/2;

        
        // Create a dictionary for this actuatorLineElement
        dictionary dict;
        dict.add("position", position);
        dict.add("coefficientData", coefficientData_);
        dict.add("chordLength", chordLength);
        dict.add("chordDirection", chordDirection);
        dict.add("spanLength", spanLength);
        dict.add("spanDirection", spanDirection);
        dict.add("freeStreamDirection", freeStreamDirection_);
        
        if (debug)
        {
            Info<< "Creating actuatorLineElement: " << name << endl;
            Info<< "Geometry segment index: " << geometrySegmentIndex << endl;
            Info<< "Position: " << position << endl;
            Info<< "Chord length: " << chordLength << endl;
            Info<< "Chord direction (before pitching): " << chordDirection 
                << endl;
            Info<< "Pitch (degrees): " << pitch << endl;
            Info<< "Span length: " << spanLength << endl;
            Info<< "Span direction: " << spanDirection << endl;
        }
        
        actuatorLineElement* element = new actuatorLineElement(name, dict, mesh_);
        elements_.set(i, element);
        pitch = pitch/180.0*Foam::constant::mathematical::pi;
        elements_[i].pitch(pitch);
        elements_[i].setVelocity(initialVelocity);
    }
}


void Foam::fv::actuatorLineSource::rotate
(
    vector rotationPoint, 
    vector axis,
    scalar radians
)
{
    forAll(elements_, i)
    {
        elements_[i].rotate(rotationPoint, axis, radians, true);
    }
}


void Foam::fv::actuatorLineSource::pitch(scalar radians)
{
    forAll(elements_, i)
    {
        elements_[i].pitch(radians);
    }
}


void Foam::fv::actuatorLineSource::translate(vector translationVector)
{
    forAll(elements_, i)
    {
        elements_[i].translate(translationVector);
    }
}


Foam::vector& Foam::fv::actuatorLineSource::force()
{
    return force_;
}


Foam::volVectorField& Foam::fv::actuatorLineSource::forceField()
{
    return forceField_;
}


Foam::vector Foam::fv::actuatorLineSource::moment(vector point)
{
    vector moment(0, 0, 0);
    forAll(elements_, i)
    {
        moment += elements_[i].moment(point);
    }
    return moment;
}


void Foam::fv::actuatorLineSource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    // Zero out force field
    forceField_ *= 0;

    // Read the reference density for incompressible flow
    //coeffs_.lookup("rhoRef") >> rhoRef_;
    
    // Zero the total force vector
    force_ = vector::zero;
    
    forAll(elements_, i)
    {
        elements_[i].addSup(eqn, forceField_);
        force_ += elements_[i].force();
    }
    
    Info<< "Force (per unit density) on " << name_ << ": "
        << endl << force_ << endl << endl;

    // Add source to eqn
    eqn += forceField_;
}


void Foam::fv::actuatorLineSource::addSup
(
    fvMatrix<scalar>& eqn,
    const label fieldI
)
{
    Info<< endl << "Adding turbulence from " << name_ << endl << endl;
    forAll(elements_, i)
    {
        elements_[i].addTurbulence(eqn);
    }
}


void Foam::fv::actuatorLineSource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    // Zero force field
    forceField_ *= 0;

    // Add source to eqn
    eqn += forceField_;
}

// ************************************************************************* //
