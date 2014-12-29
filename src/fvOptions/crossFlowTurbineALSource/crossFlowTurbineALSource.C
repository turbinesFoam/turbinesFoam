/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

#include "crossFlowTurbineALSource.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMatrices.H"
#include "geometricOneField.H"
#include "syncTools.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(crossFlowTurbineALSource, 0);
    addToRunTimeSelectionTable
    (
        option,
        crossFlowTurbineALSource,
        dictionary
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fv::crossFlowTurbineALSource::setFaceArea(vector& axis, const bool correct)
{
    area_ = 0.0;

    static const scalar tol = 0.8;

    const label nInternalFaces = mesh_.nInternalFaces();
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
    const vectorField& Sf = mesh_.Sf();
    const scalarField& magSf = mesh_.magSf();

    vector n = vector::zero;

    // calculate cell addressing for selected cells
    labelList cellAddr(mesh_.nCells(), -1);
    UIndirectList<label>(cellAddr, cells_) = identity(cells_.size());
    labelList nbrFaceCellAddr(mesh_.nFaces() - nInternalFaces, -1);
    forAll(pbm, patchI)
    {
        const polyPatch& pp = pbm[patchI];

        if (pp.coupled())
        {
            forAll(pp, i)
            {
                label faceI = pp.start() + i;
                label nbrFaceI = faceI - nInternalFaces;
                label own = mesh_.faceOwner()[faceI];
                nbrFaceCellAddr[nbrFaceI] = cellAddr[own];
            }
        }
    }

    // correct for parallel running
    syncTools::swapBoundaryFaceList(mesh_, nbrFaceCellAddr);

    // add internal field contributions
    for (label faceI = 0; faceI < nInternalFaces; faceI++)
    {
        const label own = cellAddr[mesh_.faceOwner()[faceI]];
        const label nbr = cellAddr[mesh_.faceNeighbour()[faceI]];

        if ((own != -1) && (nbr == -1))
        {
            vector nf = Sf[faceI]/magSf[faceI];

            if ((nf & axis) > tol)
            {
                area_[own] += magSf[faceI];
                n += Sf[faceI];
            }
        }
        else if ((own == -1) && (nbr != -1))
        {
            vector nf = Sf[faceI]/magSf[faceI];

            if ((-nf & axis) > tol)
            {
                area_[nbr] += magSf[faceI];
                n -= Sf[faceI];
            }
        }
    }


    // add boundary contributions
    forAll(pbm, patchI)
    {
        const polyPatch& pp = pbm[patchI];
        const vectorField& Sfp = mesh_.Sf().boundaryField()[patchI];
        const scalarField& magSfp = mesh_.magSf().boundaryField()[patchI];

        if (pp.coupled())
        {
            forAll(pp, j)
            {
                const label faceI = pp.start() + j;
                const label own = cellAddr[mesh_.faceOwner()[faceI]];
                const label nbr = nbrFaceCellAddr[faceI - nInternalFaces];
                const vector nf = Sfp[j]/magSfp[j];

                if ((own != -1) && (nbr == -1) && ((nf & axis) > tol))
                {
                    area_[own] += magSfp[j];
                    n += Sfp[j];
                }
            }
        }
        else
        {
            forAll(pp, j)
            {
                const label faceI = pp.start() + j;
                const label own = cellAddr[mesh_.faceOwner()[faceI]];
                const vector nf = Sfp[j]/magSfp[j];

                if ((own != -1) && ((nf & axis) > tol))
                {
                    area_[own] += magSfp[j];
                    n += Sfp[j];
                }
            }
        }
    }

    if (correct)
    {
        reduce(n, sumOp<vector>());
        axis = n/mag(n);
    }

    if (debug)
    {
        volScalarField area
        (
            IOobject
            (
                name_ + ":area",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("0", dimArea, 0)
        );
        UIndirectList<scalar>(area.internalField(), cells_) = area_;

        Info<< type() << ": " << name_ << " writing field " << area.name()
            << endl;

        area.write();
    }
}


void Foam::fv::crossFlowTurbineALSource::createCoordinateSystem()
{
    // construct the local rotor coordinate system
    vector origin(vector::zero);
    vector axis(vector::zero);
    vector refDir(vector::zero);

    coeffs_.lookup("origin") >> origin;
    coeffs_.lookup("axis") >> axis;

    localAxesRotation_.reset
    (
        new localAxesRotation
        (
            mesh_,
            axis,
            origin,
            cells_
        )
    );

    setFaceArea(axis, false);

    coordSys_ = cylindricalCS("rotorCoordSys", origin, axis, refDir, false);

    const scalar sumArea = gSum(area_);
    const scalar diameter = rotorRadius_*2;
    Info<< "    Rotor gometry:" << nl
        << "    - rotor diameter = " << diameter << nl
        << "    - frontal area   = " << sumArea << nl
        << "    - origin         = " << coordSys_.origin() << nl
        << "    - r-axis         = " << coordSys_.R().e1() << nl
        << "    - theta-axis     = " << coordSys_.R().e2() << nl
        << "    - z-axis         = " << coordSys_.R().e3() << endl;
}


void Foam::fv::crossFlowTurbineALSource::constructGeometry()
{
    const vectorField& C = mesh_.C();

    forAll(cells_, i)
    {
        if (area_[i] > ROOTVSMALL)
        {
            const label cellI = cells_[i];

            // position in (planar) rotor co-ordinate system
            x_[i] = coordSys_.localPosition(C[cellI]);

            // cache max radius
            rMax_ = max(rMax_, x_[i].x());

            // Wrong rotation tensor
            R_[i] = tensor(0, 0, 0, 0, 1, 0, 0, 0, 0);
            invR_[i] = R_[i].T();
        }
    }
}


Foam::tmp<Foam::vectorField> Foam::fv::crossFlowTurbineALSource::inflowVelocity
(
    const volVectorField& U
) const
{
    return U.internalField();
}


void Foam::fv::crossFlowTurbineALSource::createBlades()
{
    int nBlades = nBlades_;
    blades_.setSize(nBlades);
    int nElements;
    word profileName;
    List<List<scalar> > elementData;
    List<List<scalar> > profileData;
    word modelType = "actuatorLineSource";
    vector freeStreamDirection = freeStreamVelocity_/mag(freeStreamVelocity_);
    
    const dictionary& profilesSubDict(coeffs_.subDict("profiles"));
    
    for (int i = 0; i < nBlades_; i++)
    {
        word& bladeName = bladeNames_[i];
        // Create dictionary items for this blade
        dictionary bladeSubDict;
        bladeSubDict = bladesDict_.subDict(bladeName);
        bladeSubDict.lookup("nElements") >> nElements;
        bladeSubDict.lookup("profile") >> profileName;
        bladeSubDict.lookup("elementData") >> elementData;
        profilesSubDict.subDict(profileName).lookup("data") >> profileData;
        
        bladeSubDict.add("fieldNames", coeffs_.lookup("fieldNames"));
        bladeSubDict.add("coefficientData", profileData);
        bladeSubDict.add("tipEffect", tipEffect_);
        bladeSubDict.add("freeStreamVelocity", freeStreamVelocity_);
        
        if (debug)
        {
            Info<< "Creating actuator line blade " << bladeName << endl;
            Info<< "Blade has " << nElements << " elements" << endl;
            Info<< "Blade profile: " << profileName << endl;
            Info<< "Element data:" << endl;
            Info<< elementData << endl << endl;
            Info<< "Profile sectional coefficient data:" << endl;
            Info<< profileData << endl << endl;
        }
        
        // Convert element data into actuator line element geometry
        label nGeomPoints = elementData.size();
        List<List<List<scalar> > > elementGeometry;
        elementGeometry.setSize(nGeomPoints);
        for (int j = 0; j < nGeomPoints; j++)
        {
            // Read CFTAL dict data
            scalar axialDistance = elementData[j][0];
            scalar radius = elementData[j][1];
            scalar chordLength = elementData[j][3];
            scalar chordMount = elementData[j][4];
            
            // Set sizes for actuatorLineSource elementGeometry lists
            elementGeometry[j].setSize(4);
            elementGeometry[j][0].setSize(3);
            elementGeometry[j][1].setSize(3);
            elementGeometry[j][2].setSize(1);
            elementGeometry[j][3].setSize(1);
            
            // Create geometry point for AL source at origin
            vector point = origin_;
            // Move along axis
            point += axis_*axialDistance;
            scalar chordDisplacement = (0.5 - chordMount)*chordLength;
            point += chordDisplacement*freeStreamDirection;
            elementGeometry[j][0][0] = point.x(); // x location of geom point
            elementGeometry[j][0][1] = point.y(); // y location of geom point
            elementGeometry[j][0][2] = point.z(); // z location of geom point
            
            // Set span directions for AL source
            elementGeometry[j][1][0] = axis_.x(); // x component of span direction
            elementGeometry[j][1][1] = axis_.y(); // y component of span direction
            elementGeometry[j][1][2] = axis_.z(); // z component of span direction
            
            // Set chord length
            elementGeometry[j][2][0] = chordLength;
            
            // Set pitch
            scalar pitch = elementData[j][5];
            elementGeometry[j][3][0] = pitch;
        }
        
        if (debug)
        {
            Info<< "Converted element geometry:" << endl << elementGeometry 
                << endl;
        }
        
        bladeSubDict.add("elementGeometry", elementGeometry);
        
        dictionary dict;
        dict.add("actuatorLineSourceCoeffs", bladeSubDict);
        dict.add("type", "actuatorLineSource");
        dict.add("active", dict_.lookup("active"));
        dict.add("selectionMode", dict_.lookup("selectionMode"));
        dict.add("cellSet", dict_.lookup("cellSet"));
        
        actuatorLineSource* blade = new actuatorLineSource
        (
            bladeName, 
            modelType, 
            dict, 
            mesh_
        );
        
        blades_.set(i, blade);
    }
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::fv::crossFlowTurbineALSource::crossFlowTurbineALSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    option(name, modelType, dict, mesh),
    time_(mesh.time()),
    rhoRef_(1.0),
    omega_(0.0),
    nBlades_(0),
    freeStreamVelocity_(vector::zero),
    tipEffect_(1.0),
    x_(cells_.size(), vector::zero),
    R_(cells_.size(), I),
    invR_(cells_.size(), I),
    area_(cells_.size(), 0.0),
    coordSys_(false),
    localAxesRotation_(),
    rMax_(0.0)
{
    read(dict);
    createBlades();
    createCoordinateSystem();
    constructGeometry();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::crossFlowTurbineALSource::~crossFlowTurbineALSource()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class RhoFieldType>
void Foam::fv::crossFlowTurbineALSource::calculate
(
    const RhoFieldType& rho,
    const vectorField& U,
    vectorField& force,
    const bool divideVolume,
    const bool output
) const
{
    const scalarField& V = mesh_.V();

    // logging info
    scalar dragEff = 0.0;
    scalar liftEff = 0.0;
    scalar AOAmin = GREAT;
    scalar AOAmax = -GREAT;

    forAll(cells_, i)
    {
        if (area_[i] > ROOTVSMALL)
        {
            const label cellI = cells_[i];

            const scalar radius = x_[i].x();

            // velocity in local cylindrical reference frame
            vector Uc = localAxesRotation_->transform(U[cellI], i);

            // transform from rotor cylindrical into local coning system
            Uc = R_[i] & Uc;

            // set radial component of velocity to zero
            Uc.x() = 0.0;

            // set blade normal component of velocity
            Uc.y() = radius*omega_ - Uc.y();

            // determine blade data for this radius
            // i2 = index of upper radius bound data point in blade list
            scalar twist = 0.0;
            scalar chord = 0.0;
            scalar invDr = 0.0;

            // flip geometric angle if blade is spinning in reverse (clockwise)
            scalar alphaGeom = twist;
            if (omega_ < 0)
            {
                alphaGeom = mathematical::pi - alphaGeom;
            }

            // effective angle of attack
            scalar alphaEff = alphaGeom - atan2(-Uc.z(), Uc.y());
            if (alphaEff > mathematical::pi)
            {
                alphaEff -= mathematical::twoPi;
            }
            if (alphaEff < -mathematical::pi)
            {
                alphaEff += mathematical::twoPi;
            }

            AOAmin = min(AOAmin, alphaEff);
            AOAmax = max(AOAmax, alphaEff);

            scalar Cd1 = 0.0;
            scalar Cl1 = 0.0;

            scalar Cd2 = 0.0;
            scalar Cl2 = 0.0;

            scalar Cd = invDr*(Cd2 - Cd1) + Cd1;
            scalar Cl = invDr*(Cl2 - Cl1) + Cl1;

            // apply tip effect for blade lift
            scalar tipFactor = neg(radius/rMax_ - tipEffect_);

            // calculate forces perpendicular to blade
            scalar pDyn = 0.5*rho[cellI]*magSqr(Uc);

            scalar f = pDyn*chord*nBlades_*area_[i]/radius/mathematical::twoPi;
            vector localForce = vector(0.0, -f*Cd, tipFactor*f*Cl);

            // accumulate forces
            dragEff += rhoRef_*localForce.y();
            liftEff += rhoRef_*localForce.z();

            // convert force from local coning system into rotor cylindrical
            localForce = invR_[i] & localForce;

            // convert force to global cartesian co-ordinate system
            force[cellI] = localAxesRotation_->invTransform(localForce, i);

            if (divideVolume)
            {
                force[cellI] /= V[cellI];
            }
        }
    }

    if (output)
    {
        reduce(AOAmin, minOp<scalar>());
        reduce(AOAmax, maxOp<scalar>());
        reduce(dragEff, sumOp<scalar>());
        reduce(liftEff, sumOp<scalar>());

        Info<< type() << " output:" << nl
            << "    min/max(AOA)   = " << radToDeg(AOAmin) << ", "
            << radToDeg(AOAmax) << nl
            << "    Effective drag = " << dragEff << nl
            << "    Effective lift = " << liftEff << endl;
    }
}


void Foam::fv::crossFlowTurbineALSource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    volVectorField force
    (
        IOobject
        (
            name_ + ":rotorForce",
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
    coeffs_.lookup("rhoRef") >> rhoRef_;

    const vectorField Uin(inflowVelocity(eqn.psi()));
    calculate(geometricOneField(), Uin, force);

    // Add source to rhs of eqn
    eqn -= force;

    if (mesh_.time().outputTime())
    {
        force.write();
    }
}


void Foam::fv::crossFlowTurbineALSource::addSup
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
            name_ + ":rotorForce",
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

    const vectorField Uin(inflowVelocity(eqn.psi()));
    calculate(rho, Uin, force);

    // Add source to rhs of eqn
    eqn -= force;

    if (mesh_.time().outputTime())
    {
        force.write();
    }
}


void Foam::fv::crossFlowTurbineALSource::writeData(Ostream& os) const
{
    os  << indent << name_ << endl;
    dict_.write(os);
}


void Foam::fv::crossFlowTurbineALSource::printCoeffs() const
{
    Info<< "Number of blades: " << nBlades_ << endl;
}


bool Foam::fv::crossFlowTurbineALSource::read(const dictionary& dict)
{
    if (option::read(dict))
    {
        
        coeffs_.lookup("fieldNames") >> fieldNames_;
        applied_.setSize(fieldNames_.size(), false);

        // Read coordinate system/geometry invariant properties
        coeffs_.lookup("origin") >> origin_;
        coeffs_.lookup("axis") >> axis_;
        coeffs_.lookup("freeStreamVelocity") >> freeStreamVelocity_;
        coeffs_.lookup("tipSpeedRatio") >> tipSpeedRatio_;
        coeffs_.lookup("rotorRadius") >> rotorRadius_;
        omega_ = tipSpeedRatio_*mag(freeStreamVelocity_)/rotorRadius_;

        // Get blade information
        bladesDict_ = coeffs_.subDict("blades");
        nBlades_ = bladesDict_.keys().size();
        bladeNames_ = bladesDict_.toc();

        coeffs_.lookup("tipEffect") >> tipEffect_;
        
        if (debug)
        {
            Info<< "Debugging on" << endl;
            Info<< "Cross-flow turbine properties:" << endl;
            printCoeffs();
        }

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
