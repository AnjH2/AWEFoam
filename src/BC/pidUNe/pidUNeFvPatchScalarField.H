/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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


See also
    Foam::fixedValueFvPatchField

SourceFiles
    pidUNeFvPatchScalarField.C

======================== pidUNeFvPatchScalarField.H ========================*/
#ifndef pidUNeFvPatchScalarField_H
#define pidUNeFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"
#include "Enum.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "fvMesh.H"         // usually already indirectly included, but safe


namespace Foam
{

class pidUNeFvPatchScalarField : public fixedValueFvPatchScalarField
{
public:
    enum controlMode { cmTotal, cmAverage };
    static const Enum<controlMode> controlModeNames_;

private:
    // Controller params
    controlMode mode_;
    word   currentFieldName_;
    scalar jGoal_;
    scalar Kp_, Ki_;
    scalar dUmax_;
    scalar Umin_, Umax_;
    label  updateEveryNTimeSteps_;

    // State
    label  lastUpdateTimeIndex_;
    bool   initFromPatchAverage_;
    bool   initialized_;
    scalar U_;
    scalar eInt_;

    // Printing
    bool   verbose_;
    label  lastPrintTimeIndex_;

    // External dict support (system/pidUNeDict)
    bool                 useExternalDict_;
    word                 externalDictName_;   // e.g. "pidUNeDict"
    autoPtr<IOdictionary> extDictPtr_;
    label extDictEventNo_;

    // Helpers
    void ensureExternalDict();
    void applyDict(const dictionary& d);
    const dictionary* externalDictTop() const;
    void readAndApplyExternalDict();

public:
    TypeName("pidUNe");

    pidUNeFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    pidUNeFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const dictionary&
    );

    pidUNeFvPatchScalarField(const pidUNeFvPatchScalarField&);
    pidUNeFvPatchScalarField
    (
    const pidUNeFvPatchScalarField&,
    const DimensionedField<scalar, volMesh>&
);

    pidUNeFvPatchScalarField
    (
        const pidUNeFvPatchScalarField&,
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const fvPatchFieldMapper&
    );

    virtual tmp<fvPatchScalarField> clone() const
    {
        return tmp<fvPatchScalarField>(new pidUNeFvPatchScalarField(*this));
    }

    virtual tmp<fvPatchScalarField> clone
    (
        const DimensionedField<scalar, volMesh>& iF
    ) const
    {
        return tmp<fvPatchScalarField>(new pidUNeFvPatchScalarField(*this, iF));
    }

    virtual void updateCoeffs();
    virtual void write(Ostream& os) const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
