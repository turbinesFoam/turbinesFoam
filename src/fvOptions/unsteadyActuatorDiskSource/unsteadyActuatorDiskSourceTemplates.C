/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020 ENERCON GmbH
    Copyright (C) 2018-2020 OpenCFD Ltd
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

#include "unsteadyActuatorDiskSource.H"
#include "fvMesh.H"
#include "fvMatrix.H"
#include "volFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class AlphaFieldType, class RhoFieldType>
void Foam::fv::unsteadyActuatorDiskSource::calc
(
    const AlphaFieldType& alpha,
    const RhoFieldType& rho,
    fvMatrix<vector>& eqn
)
{
    switch (forceMethod_)
    {
        case forceMethodType::FROUDE:
        {
            calcFroudeMethod(alpha, rho, eqn);
            break;
        }

        case forceMethodType::VARIABLE_SCALING:
        {
            calcVariableScalingMethod(alpha, rho, eqn);
            break;
        }

        default:
            break;
    }
}


template<class AlphaFieldType, class RhoFieldType>
void Foam::fv::unsteadyActuatorDiskSource::calcFroudeMethod
(
    const AlphaFieldType& alpha,
    const RhoFieldType& rho,
    fvMatrix<vector>& eqn
)
{
    const vectorField& U = eqn.psi();
    vectorField& Usource = eqn.source();
    const scalarField& cellsV = mesh_.V();

    // Compute upstream U and rho, spatial-averaged over monitor-region
    vector Uref(Zero);
    scalar rhoRef = 0.0;
    label szMonitorCells = monitorCells_.size();

    for (const label celli : monitorCells_)
    {
        Uref += U[celli];
        rhoRef = rhoRef + rho[celli];
    }
    reduce(Uref, sumOp<vector>());
    reduce(rhoRef, sumOp<scalar>());
    reduce(szMonitorCells, sumOp<label>());

    if (szMonitorCells == 0)
    {
        FatalErrorInFunction
            << "No cell is available for incoming velocity monitoring."
            << exit(FatalError);
    }

    Uref /= szMonitorCells;
    rhoRef /= szMonitorCells;

    const scalar Ct = sink_*UvsCtPtr_->value(mag(Uref));
    const scalar Cp = sink_*UvsCpPtr_->value(mag(Uref));

    if (Cp <= VSMALL || Ct <= VSMALL)
    {
        FatalErrorInFunction
           << "Cp and Ct must be greater than zero." << nl
           << "Cp = " << Cp << ", Ct = " << Ct
           << exit(FatalError);
    }

    // (BJSB:Eq. 3.9)
    const vector diskDir = this->diskDir();
    const scalar a = 1.0 - Cp/Ct;
    const scalar T = 2.0*rhoRef*diskArea_*magSqr(Uref & diskDir)*a*(1 - a);

    for (const label celli : cells_)
    {
        Usource[celli] += ((cellsV[celli]/V())*T)*diskDir;
    }

    if
    (
        mesh_.time().timeOutputValue() >= writeFileStart_
     && mesh_.time().timeOutputValue() <= writeFileEnd_
    )
    {
        Ostream& os = file();
        writeCurrentTime(os);

        os  << Uref << tab << Cp << tab << Ct << tab << a << tab << T
            << endl;
    }
}


template<class AlphaFieldType, class RhoFieldType>
void Foam::fv::unsteadyActuatorDiskSource::calcVariableScalingMethod
(
    const AlphaFieldType& alpha,
    const RhoFieldType& rho,
    fvMatrix<vector>& eqn
)
{
    const vectorField& U = eqn.psi();
    vectorField& Usource = eqn.source();
    const scalarField& cellsV = mesh_.V();

    // Monitor and average monitor-region U and rho
    vector Uref(Zero);
    scalar rhoRef = 0.0;
    label szMonitorCells = monitorCells_.size();

    for (const label celli : monitorCells_)
    {
        Uref += U[celli];
        rhoRef = rhoRef + rho[celli];
    }
    reduce(Uref, sumOp<vector>());
    reduce(rhoRef, sumOp<scalar>());
    reduce(szMonitorCells, sumOp<label>());

    if (szMonitorCells == 0)
    {
        FatalErrorInFunction
            << "No cell is available for incoming velocity monitoring."
            << exit(FatalError);
    }

    Uref /= szMonitorCells;
    const scalar magUref = mag(Uref);
    rhoRef /= szMonitorCells;

    // Monitor and average U and rho on actuator disk
    vector Udisk(Zero);
    scalar rhoDisk = 0.0;
    scalar totalV = 0.0;

    for (const label celli : cells_)
    {
        Udisk += U[celli]*cellsV[celli];
        rhoDisk += rho[celli]*cellsV[celli];
        totalV += cellsV[celli];
    }
    reduce(Udisk, sumOp<vector>());
    reduce(rhoDisk, sumOp<scalar>());
    reduce(totalV, sumOp<scalar>());

    if (totalV < SMALL)
    {
        FatalErrorInFunction
            << "No cell in the actuator disk."
            << exit(FatalError);
    }

    Udisk /= totalV;
    const scalar magUdisk = mag(Udisk);
    rhoDisk /= totalV;

    if (mag(Udisk) < SMALL)
    {
        FatalErrorInFunction
            << "Velocity spatial-averaged on actuator disk is zero." << nl
            << "Please check if the initial U field is zero."
            << exit(FatalError);
    }

    // Interpolated thrust/power coeffs from power/thrust curves
    const scalar Ct = sink_*UvsCtPtr_->value(magUref);
    const scalar Cp = sink_*UvsCpPtr_->value(magUref);

    if (Cp <= VSMALL || Ct <= VSMALL)
    {
        FatalErrorInFunction
           << "Cp and Ct must be greater than zero." << nl
           << "Cp = " << Cp << ", Ct = " << Ct
           << exit(FatalError);
    }

    // Calibrated thrust/power coeffs from power/thrust curves (LSRMTK:Eq. 6)
    const scalar CtStar = Ct*sqr(magUref/magUdisk);
    const scalar CpStar = Cp*pow3(magUref/magUdisk);

    // Compute calibrated thrust/power (LSRMTK:Eq. 5)
    const vector diskDir = this->diskDir();
    const scalar T = 0.5*rhoRef*diskArea_*magSqr(Udisk & diskDir)*CtStar;
    const scalar P = 0.5*rhoRef*diskArea_*pow3(mag(Udisk & diskDir))*CpStar;

    for (const label celli : cells_)
    {
        Usource[celli] += (cellsV[celli]/totalV*T)*diskDir;
    }

    if
    (
        mesh_.time().timeOutputValue() >= writeFileStart_
     && mesh_.time().timeOutputValue() <= writeFileEnd_
    )
    {
        Ostream& os = file();
        writeCurrentTime(os);

        os  << Uref << tab << Cp << tab << Ct << tab
            << Udisk << tab << CpStar << tab << CtStar << tab << T << tab << P
            << endl;
    }
}


// ************************************************************************* //
