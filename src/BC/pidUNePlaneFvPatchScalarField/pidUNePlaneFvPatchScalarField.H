/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

Class
    Foam::pidUNePlaneFvPatchScalarField

Description
    Fixed-value scalar boundary condition for UNe with a PI controller.

    This is the plane/faceZone variant of pidUNeFvPatchScalarField.  The
    measured current is calculated through a named faceZone, using the same
    projected-plane integration style as pidUNeConstraint:

        I = integral (J dot nHat) dA_projected

    when useProjected is true, otherwise the actual face-normal flux
    integral integral (J dot Sf) is used.

SourceFiles
    pidUNePlaneFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef pidUNePlaneFvPatchScalarField_H
#define pidUNePlaneFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"
#include "Enum.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "fvMesh.H"

namespace Foam
{

class pidUNePlaneFvPatchScalarField : public fixedValueFvPatchScalarField
{
public:
    enum controlMode { cmTotal, cmAverage };
    static const Enum<controlMode> controlModeNames_;

private:
    // Controller params
    controlMode mode_;
    word   currentFieldName_;
    word   faceZoneName_;
    bool   useProjected_;
    vector nHat_;

    scalar jGoal_;
    scalar Kp_, Ki_;
    scalar dUmax_;
    scalar Umin_, Umax_;
    label  updateEveryNTimeSteps_;
    bool   firstControllerUpdateSkipped_;

    // State
    label  lastUpdateTimeIndex_;
    bool   initFromPatchAverage_;
    bool   initialized_;
    scalar U_;
    scalar eInt_;

    // Printing
    bool   verbose_;
    label  lastPrintTimeIndex_;

    // External dict support (system/pidUNeDict by default)
    bool                  useExternalDict_;
    word                  externalDictName_;
    autoPtr<IOdictionary> extDictPtr_;
    label                 extDictEventNo_;

    // Helpers
    void ensureExternalDict();
    void applyDict(const dictionary& d);
    const dictionary* externalDictTop() const;
    void readAndApplyExternalDict();

    void calculatePlaneCurrent
    (
        scalar& Iabs,
        scalar& A,
        scalar& jAvg
    ) const;

public:
    TypeName("pidUNePlane");

    pidUNePlaneFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    pidUNePlaneFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const dictionary&
    );

    pidUNePlaneFvPatchScalarField(const pidUNePlaneFvPatchScalarField&);

    pidUNePlaneFvPatchScalarField
    (
        const pidUNePlaneFvPatchScalarField&,
        const DimensionedField<scalar, volMesh>&
    );

    pidUNePlaneFvPatchScalarField
    (
        const pidUNePlaneFvPatchScalarField&,
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const fvPatchFieldMapper&
    );

    virtual tmp<fvPatchScalarField> clone() const
    {
        return tmp<fvPatchScalarField>
        (
            new pidUNePlaneFvPatchScalarField(*this)
        );
    }

    virtual tmp<fvPatchScalarField> clone
    (
        const DimensionedField<scalar, volMesh>& iF
    ) const
    {
        return tmp<fvPatchScalarField>
        (
            new pidUNePlaneFvPatchScalarField(*this, iF)
        );
    }

    virtual void updateCoeffs();
    virtual void write(Ostream& os) const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
