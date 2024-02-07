/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020 ENERCON GmbH
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::fv::unsteadyActuatorDiskSource

Group
    grpFvOptionsSources

Description
    Applies sources on velocity (i.e. \c U) within a specified region
    to enable actuator disk models for aero/hydro thrust loading of
    horizontal axis turbines on surrounding flow field in terms of
    energy conversion processes.

    Available options for force computations:
    \verbatim
      Froude          | Froude's one-dimensional ideal actuator disk method
      variableScaling | Variable-scaling actuator disk method
    \endverbatim

    The expressions for \c Froude method:
    \f[
        T = 2 \rho A | \vec{U_o} \cdot \vec{n} |^2 a (1-a)
    \f]
    with
    \f[
        a = 1 - \frac{C_p}{C_T}
    \f]
    where
    \vartable
      T         | Thrust magnitude
      \rho      | Monitored incoming fluid density
      A         | Actuator disk planar surface area
      \vec{U_o} | Incoming velocity spatial-averaged on monitored region
      \vec{n}   | Surface-normal vector of the actuator disk pointing downstream
      a         | Axial induction factor
      C_p       | Power coefficient
      C_T       | Thrust coefficient
    \endvartable

    The expressions for \c variableScaling method:
    \f[
        T = 0.5 \rho A | \vec{U_m} \cdot \vec{n} |^2 C_T^*
    \f]
    with
    \f[
        C_T^* = C_T \left( \frac{ |\vec{U_{ref}}| }{ |\vec{U_m}| } \right)^2
    \f]
    where
    \vartable
      T         | Thrust magnitude
      \rho      | Monitored incoming fluid density
      A         | Actuator disk planar surface area
      \vec{U_m} | Incoming velocity spatial-averaged on actuator disk
      \vec{n}   | Surface-normal vector of the actuator disk pointing downstream
      \vec{U_{ref}} | Incoming velocity spatial-averaged on monitored region
      C_T       | Thrust coefficient
      C_T^*     | Calibrated thrust coefficient
    \endvartable

    Reference
    \verbatim
        Froude's one-dimensional ideal actuator disk (tags:F,BJSB):
            Froude, R. E. (1889).
            On the part played in propulsion by differences of fluid pressure.
            Trans. Inst. Naval Architects, 30, 390.

            Burton, T., Jenkins, N., Sharpe, D., & Bossanyi, E. (2011).
            Wind energy handbook, 2nd edition.
            Chichester, the United Kingdom. John Wiley & Sons.
            Print ISBN:978-0-470-69975-1
            DOI:10.1002/9781119992714

        Variable-scaling method (tags:LSRMKT,LS):
            van der Laan, M. P., Sørensen, N. N., Réthoré, P. E.,
            Mann, J., Kelly, M. C., & Troldborg, N. (2015).
            The k‐ε‐fP model applied to double wind turbine
            wakes using different actuator disk force methods.
            Wind Energy, 18(12), 2223-2240.
            DOI:10.1002/we.1816

            van der Laan, M. P., Sørensen, N. N., Réthoré, P. E., Mann, J.,
            Kelly, M. C., Troldborg, N., Hansen, K. S., & Murcia, J. P. (2015).
            The k‐ε‐fP model applied to wind farms.
            Wind Energy, 18(12), 2065-2084.
            DOI:10.1002/we.1804
    \endverbatim

    Sources applied to:
    \verbatim
      U         | Velocity    [m/s]
    \endverbatim

    Required fields:
    \verbatim
      U         | Velocity    [m/s]
    \endverbatim

Usage
    Minimal example by using \c constant/fvOptions:
    \verbatim
    unsteadyActuatorDiskSource1
    {
        // Mandatory entries (unmodifiable)
        type            unsteadyActuatorDiskSource;

        // Mandatory entries (runtime modifiable)
        diskArea        40.0;
        diskDir         (1 0 0);
        Cp              <Function1>;
        Ct              <Function1>;

        // Conditional optional entries (unmodifiable)
        monitorMethod   <points|cellSet>;
        monitorCoeffs
        {
            // Option-1
            points
            (
                (p1x p1y p1z)
                (p2x p2y p2z)
                ...
            );

            // Option-2
            cellSet     <monitorCellSet>;
        }

        // Optional entries (unmodifiable)
        variant         <forceMethod>;

        // Optional entries (runtime modifiable)
        sink            true;
        writeFileStart  0;
        writeFileEnd    100;

        // Mandatory/Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property   | Description                         | Type | Reqd     | Dflt
      type       | Type name: unsteadyActuatorDiskSource      | word | yes      | -
      diskArea   | Actuator disk planar surface area   | scalar | yes    | -
      diskDir    | Surface-normal vector of the actuator disk <!--
               --> pointing downstream                 | vector    | yes | -
      Cp         | Power coefficient                   | Function1 | yes | -
      Ct         | Thrust coefficient                  | Function1 | yes | -
      monitorMethod | Type of incoming velocity monitoring method   <!--
               -->    - see below                      | word | no  | points
      variant    | Type of the force computation method - see below <!--
               -->                                     | word | no  | Froude
      sink       | Flag for body forces to act as a source (false)  <!--
               -->  or a sink (true)                   | bool   | no  | true
      writeFileStart | Start time for file output      | scalar | no  | 0
      writeFileEnd   | End time for file output        | scalar | no  | VGREAT
    \endtable

    The inherited entries are elaborated in:
      - \link fvOption.H \endlink
      - \link cellSetOption.H \endlink
      - \link writeFile.H \endlink
      - \link Function1.H \endlink

    Options for the \c monitorMethod entry:
    \verbatim
      points  | Monitor incoming velocity field at a given set of points
      cellSet | Monitor incoming velocity field at a given cellSet
    \endverbatim

    Options for the \c variant entry:
    \verbatim
      Froude          | Froude's one-dimensional ideal actuator disk method
      variableScaling | Variable-scaling actuator disk method
    \endverbatim

SourceFiles
    unsteadyActuatorDiskSource.C
    unsteadyActuatorDiskSourceTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef unsteadyActuatorDiskSource_H
#define unsteadyActuatorDiskSource_H

#include "cellSetOption.H"
#include "writeFile.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                     Class unsteadyActuatorDiskSource Declaration
\*---------------------------------------------------------------------------*/

class unsteadyActuatorDiskSource
:
    public fv::cellSetOption,
    public functionObjects::writeFile
{
protected:

    // Protected Enumerations

        //- Options for the force computation method types
        enum forceMethodType
        {
            FROUDE,                 //!< "Froude's ideal actuator disk method"
            VARIABLE_SCALING,       //!< "Variable-scaling actuator disk method"
        };

        //- Names for forceMethodType
        static const Enum<forceMethodType> forceMethodTypeNames;

        //- Options for the incoming velocity monitoring method types
        enum monitorMethodType
        {
            POINTS,                 //!< "A set of points"
            CELLSET                 //!< "A cellSet"
        };

        //- Names for monitorMethodType
        static const Enum<monitorMethodType> monitorMethodTypeNames;


    // Protected Data

        //- The type of the force computation method
        const enum forceMethodType forceMethod_;

        //- The type of incoming velocity monitoring method
        const enum monitorMethodType monitorMethod_;

        //- Flag for body forces to act as a source (false) or a sink (true)
        label sink_;

        //- Start time for file output
        scalar writeFileStart_;

        //- End time for file output
        scalar writeFileEnd_;

        //- Actuator disk planar surface area [m2]
        scalar diskArea_;

        //- Surface-normal vector of the actuator disk pointing downstream
        autoPtr<Function1<vector>> diskDir_;

        //- Velocity vs power coefficients
        autoPtr<Function1<scalar>> UvsCpPtr_;

        //- Velocity vs thrust coefficients
        autoPtr<Function1<scalar>> UvsCtPtr_;

        //- Set of cells whereat the incoming velocity is monitored
        labelList monitorCells_;


    // Protected Member Functions

        //- Output file header information
        virtual void writeFileHeader(Ostream& os);


private:

    // Private Member Functions

        //- Locate the set of cells whereat the incoming velocity is monitored
        void setMonitorCells(const dictionary& dict);

        //- Compute the selected force computation method type
        template<class AlphaFieldType, class RhoFieldType>
        void calc
        (
            const AlphaFieldType& alpha,
            const RhoFieldType& rho,
            fvMatrix<vector>& eqn
        );

        //- Compute Froude's ideal actuator disk method
        template<class AlphaFieldType, class RhoFieldType>
        void calcFroudeMethod
        (
            const AlphaFieldType& alpha,
            const RhoFieldType& rho,
            fvMatrix<vector>& eqn
        );

        //- Compute variable-scaling actuator disk method
        template<class AlphaFieldType, class RhoFieldType>
        void calcVariableScalingMethod
        (
            const AlphaFieldType& alpha,
            const RhoFieldType& rho,
            fvMatrix<vector>& eqn
        );


public:

    //- Runtime type information
    TypeName("unsteadyActuatorDiskSource");


    // Constructors

        //- No default construct
        unsteadyActuatorDiskSource() = delete;

        //- Construct from components
        unsteadyActuatorDiskSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- No copy construct
        unsteadyActuatorDiskSource(const unsteadyActuatorDiskSource&) = delete;

        //- No copy assignment
        void operator=(const unsteadyActuatorDiskSource&) = delete;


    //- Destructor
    virtual ~unsteadyActuatorDiskSource() = default;


    // Member Functions

        // Access

            //- Normal disk direction
            vector diskDir() const;

            //- Disk area
            scalar diskArea() const
            {
                return diskArea_;
            }


         // Add explicit and implicit contributions

            //- Source term to momentum equation
            virtual void addSup
            (
                fvMatrix<vector>& eqn,
                const label fieldi
            );

            //- Source term to compressible momentum equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                const label fieldi
            );

            //- Source term to phase-compressible momentum equation
            virtual void addSup
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                const label fieldi
            );


        // IO

            //- Read dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "unsteadyActuatorDiskSourceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
