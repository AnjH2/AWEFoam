/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2017 OpenFOAM Foundation
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
    Foam::fv::myFixedTemperatureConstraint

Group
    grpFvOptionsConstraints

Description
    Constrain temperature equation (i.e. \c T) with
    a given set of fixed values within a specified region.

    Constraints applied to:
    \verbatim
      T    | Temperature                   [K]
    \endverbatim

    Required fields:
    \verbatim
      T    | Temperature                   [K]
      e/h  | Internal energy/Enthalphy     [m2/s2]
    \endverbatim

Usage
    Minimal example by using \c constant/fvOptions:
    \verbatim
    myFixedTemperatureConstraint1
    {
        // Mandatory entries (unmodifiable)
        type            myFixedTemperatureConstraint;
        mode            <mode>;

        // Conditional mandatory/optional entries (runtime modifiable)

            // when mode=uniform (mandatory)
            temperature     constant 500; // fixed temperature with time [K]

            // when mode=lookup (optional)
            T            <Tname>;

        // Mandatory/Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                           | Type   | Reqd  | Dflt
      type      | Type name: myFixedTemperatureConstraint | word   | yes   | -
      mode      | Temperature mode - see below          | word   | yes   | -
      temperature | Function1 type temperature field    | Function1 | cndtnl | -
      T         | Name of operand temperature field     | word   | cndtnl | T
    \endtable

    Options for the \c mode entry:
    \verbatim
      uniform | Use Function1 type temperature field
      lookup  | Use temperature field given field name
    \endverbatim

    The inherited entries are elaborated in:
     - \link fvOption.H \endlink
     - \link cellSetOption.H \endlink

Note:
  - The \c uniform option allows the use of a time-varying
    uniform temperature by means of the \c Function1 type.

See also
  - Foam::fv::limitTemperature
  - Foam::Function1

SourceFiles
    myFixedTemperatureConstraint.C

\*---------------------------------------------------------------------------*/

#ifndef pidUNeConstraint_H
#define pidUNeConstraint_H

#include "cellSetOption.H"
#include "Enum.H"

namespace Foam
{
namespace fv
{

class pidUNeConstraint : public fv::cellSetOption
{
public:
    enum controlMode
    {
        cmTotalCurrent,
        cmAverageFlux
    };

    static const Enum<controlMode> controlModeNames_;

protected:
    controlMode mode_;

    bool   useProjected_;
    vector nHat_;                 // normalized

    word   faceZoneName_;
    word   currentFieldName_;

    scalar jGoal_, Kp_, Ki_, dUmax_, Umin_, Umax_, initialU_;
    label  Nupdate_;

    bool   verbose_;
    bool   printEveryTimeStep_;
    bool   printOnUpdate_;
    label  lastPrintTimeIndex_;

    scalar U_, eInt_;
    label  k_;
    
    
    bool initFromField_;
    bool initialized_;
    
    label updateEveryNTimeSteps_;
    label lastUpdateTimeIndex_;

public:
    TypeName("pidUNeConstraint");

    pidUNeConstraint(const word&, const word&, const dictionary&, const fvMesh&);
    virtual bool read(const dictionary& dict);
    virtual void constrain(fvMatrix<scalar>& eqn, const label fieldi);
};

} // namespace fv
} // namespace Foam

#endif


// ************************************************************************* //
