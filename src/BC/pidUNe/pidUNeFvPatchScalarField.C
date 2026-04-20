/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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


======================== pidUNeFvPatchScalarField.C ========================*/
#include "pidUNeFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvCFD.H"
#include "fvc.H"
#include "PstreamReduceOps.H"

namespace Foam
{
    makePatchTypeField(fvPatchScalarField, pidUNeFvPatchScalarField);
}

const Foam::Enum<Foam::pidUNeFvPatchScalarField::controlMode>
Foam::pidUNeFvPatchScalarField::controlModeNames_
({
    { controlMode::cmTotal,   "total"   },
    { controlMode::cmAverage, "average" }
});


void Foam::pidUNeFvPatchScalarField::applyDict(const dictionary& d)
{
    mode_ = controlModeNames_.getOrDefault("controlMode", d, mode_);

    d.readIfPresent("currentField", currentFieldName_);
    d.readIfPresent("jGoal", jGoal_);
    d.readIfPresent("Kp", Kp_);
    d.readIfPresent("Ki", Ki_);
    d.readIfPresent("dUmax", dUmax_);
    d.readIfPresent("Umin", Umin_);
    d.readIfPresent("Umax", Umax_);
    d.readIfPresent("updateEveryNTimeSteps", updateEveryNTimeSteps_);
    d.readIfPresent("initFromPatchAverage", initFromPatchAverage_);
    d.readIfPresent("verbose", verbose_);

    if (updateEveryNTimeSteps_ < 1) updateEveryNTimeSteps_ = 1;

    bool resetController(false);
    if (d.readIfPresent("resetController", resetController) && resetController)
    {
        // Optional: allow new initialU when resetting
        d.readIfPresent("initialU", U_);
        eInt_ = 0.0;
        initialized_ = false;
        lastUpdateTimeIndex_ = -1;
    }

    // Clamp immediately if limits changed
    U_ = max(Umin_, min(Umax_, U_));
}


// Return dict content (prefer registry IOdictionary if present)
const Foam::dictionary* Foam::pidUNeFvPatchScalarField::externalDictTop() const
{
    if (!useExternalDict_) return nullptr;

    const fvMesh& m = patch().boundaryMesh().mesh();

    if (m.foundObject<IOdictionary>(externalDictName_))
    {
        return &m.lookupObject<IOdictionary>(externalDictName_);
    }

    if (extDictPtr_.valid())
    {
        return &extDictPtr_();
    }

    return nullptr;
}


// Ensure we have an IOdictionary instance if none exists in registry
void Foam::pidUNeFvPatchScalarField::ensureExternalDict()
{
    if (!useExternalDict_) return;

    const fvMesh& m = patch().boundaryMesh().mesh();

    if (m.foundObject<IOdictionary>(externalDictName_)) return;
    if (extDictPtr_.valid()) return;

    extDictPtr_.reset
    (
        new IOdictionary
        (
            IOobject
            (
                externalDictName_,                 // "pidUNeDict"
                m.time().system(),                 // "system"
                m,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE
            )
        )
    );
}


void Foam::pidUNeFvPatchScalarField::readAndApplyExternalDict()
{
    if (!useExternalDict_) return;

    ensureExternalDict();

    const fvMesh& m = patch().boundaryMesh().mesh();

    // Trigger readIfModified on whichever IOdictionary we will use
    word dictPath("unknown");

    if (m.foundObject<IOdictionary>(externalDictName_))
    {
        IOdictionary& D = const_cast<IOdictionary&>
        (
            m.lookupObject<IOdictionary>(externalDictName_)
        );

        D.readIfModified();
        dictPath = D.objectPath();
    }
    else if (extDictPtr_.valid())
    {
        extDictPtr_().readIfModified();
        dictPath = extDictPtr_().objectPath();
    }

    const dictionary* topPtr = externalDictTop();
    if (!topPtr) return;

    const dictionary& top = *topPtr;

    // Select subdict: patch-name > default > top
    const dictionary* usePtr = nullptr;

    if (top.found(patch().name()))
    {
        usePtr = &top.subDict(patch().name());
    }
    else if (top.found("default"))
    {
        usePtr = &top.subDict("default");
    }
    else
    {
        usePtr = &top;
    }

    // ALWAYS apply (no eventNo logic)
    applyDict(*usePtr);

    // Optional debug
    if (verbose_ && Pstream::master())
    {
        const scalar jFile = usePtr->lookupOrDefault<scalar>("jGoal", jGoal_);
        Info<< "pidUNe: applied " << externalDictName_
            << " patch=" << patch().name()
            << " jGoal(file)=" << jFile
            << " jGoal(mem)=" << jGoal_
            << " updateN=" << updateEveryNTimeSteps_
            << " path=\"" << dictPath << "\""
            << nl;
    }
}

// Default ctor (p,iF)
Foam::pidUNeFvPatchScalarField::pidUNeFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    mode_(cmTotal),
    currentFieldName_("Ie"),
    jGoal_(0.0),
    Kp_(0.0),
    Ki_(0.0),
    dUmax_(0.1),
    Umin_(-GREAT),
    Umax_(GREAT),
    updateEveryNTimeSteps_(1),
    lastUpdateTimeIndex_(-1),
    initFromPatchAverage_(true),
    initialized_(false),
    U_(0.0),
    eInt_(0.0),
    verbose_(false),
    lastPrintTimeIndex_(-1),
    useExternalDict_(true),
    externalDictName_("pidUNeDict"),
    extDictPtr_(nullptr),
    extDictEventNo_(-1)          // IMPORTANT
{}


// Dict ctor (p,iF,dict)
Foam::pidUNeFvPatchScalarField::pidUNeFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    mode_(controlModeNames_.getOrDefault("controlMode", dict, cmTotal)),
    currentFieldName_(dict.getOrDefault<word>("currentField", "Ie")),
    jGoal_(dict.getOrDefault<scalar>("jGoal", 0.0)),
    Kp_(dict.getOrDefault<scalar>("Kp", 0.0)),
    Ki_(dict.getOrDefault<scalar>("Ki", 0.0)),
    dUmax_(dict.getOrDefault<scalar>("dUmax", 0.1)),
    Umin_(dict.getOrDefault<scalar>("Umin", -GREAT)),
    Umax_(dict.getOrDefault<scalar>("Umax", GREAT)),
    updateEveryNTimeSteps_(dict.getOrDefault<label>("updateEveryNTimeSteps", 1)),
    lastUpdateTimeIndex_(-1),
    initFromPatchAverage_(dict.getOrDefault<bool>("initFromPatchAverage", true)),
    initialized_(false),
    U_(dict.getOrDefault<scalar>("initialU", 0.0)),
    eInt_(0.0),
    verbose_(dict.getOrDefault<bool>("verbose", true)),
    lastPrintTimeIndex_(-1),
    useExternalDict_(dict.getOrDefault<bool>("useExternalDict", true)),
    externalDictName_(dict.getOrDefault<word>("externalDictName", "pidUNeDict")),
    extDictPtr_(nullptr),
    extDictEventNo_(-1)          // IMPORTANT
{
    if (updateEveryNTimeSteps_ < 1) updateEveryNTimeSteps_ = 1;

    if (!dict.found("value"))
    {
        this->operator==(U_);
    }
}


// Copy ctor
Foam::pidUNeFvPatchScalarField::pidUNeFvPatchScalarField
(
    const pidUNeFvPatchScalarField& ptf
)
:
    fixedValueFvPatchScalarField(ptf),
    mode_(ptf.mode_),
    currentFieldName_(ptf.currentFieldName_),
    jGoal_(ptf.jGoal_),
    Kp_(ptf.Kp_),
    Ki_(ptf.Ki_),
    dUmax_(ptf.dUmax_),
    Umin_(ptf.Umin_),
    Umax_(ptf.Umax_),
    updateEveryNTimeSteps_(ptf.updateEveryNTimeSteps_),
    lastUpdateTimeIndex_(ptf.lastUpdateTimeIndex_),
    initFromPatchAverage_(ptf.initFromPatchAverage_),
    initialized_(ptf.initialized_),
    U_(ptf.U_),
    eInt_(ptf.eInt_),
    verbose_(ptf.verbose_),
    lastPrintTimeIndex_(ptf.lastPrintTimeIndex_),
    useExternalDict_(ptf.useExternalDict_),
    externalDictName_(ptf.externalDictName_),
    extDictPtr_(nullptr),
    extDictEventNo_(-1)          // IMPORTANT: do NOT copy ptf.extDictEventNo_
{}


// Copy-with-iF ctor (clone(iF))
Foam::pidUNeFvPatchScalarField::pidUNeFvPatchScalarField
(
    const pidUNeFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    mode_(ptf.mode_),
    currentFieldName_(ptf.currentFieldName_),
    jGoal_(ptf.jGoal_),
    Kp_(ptf.Kp_),
    Ki_(ptf.Ki_),
    dUmax_(ptf.dUmax_),
    Umin_(ptf.Umin_),
    Umax_(ptf.Umax_),
    updateEveryNTimeSteps_(ptf.updateEveryNTimeSteps_),
    lastUpdateTimeIndex_(ptf.lastUpdateTimeIndex_),
    initFromPatchAverage_(ptf.initFromPatchAverage_),
    initialized_(ptf.initialized_),
    U_(ptf.U_),
    eInt_(ptf.eInt_),
    verbose_(ptf.verbose_),
    lastPrintTimeIndex_(ptf.lastPrintTimeIndex_),
    useExternalDict_(ptf.useExternalDict_),
    externalDictName_(ptf.externalDictName_),
    extDictPtr_(nullptr),
    extDictEventNo_(-1)          // IMPORTANT
{}


// Mapper ctor
Foam::pidUNeFvPatchScalarField::pidUNeFvPatchScalarField
(
    const pidUNeFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    mode_(ptf.mode_),
    currentFieldName_(ptf.currentFieldName_),
    jGoal_(ptf.jGoal_),
    Kp_(ptf.Kp_),
    Ki_(ptf.Ki_),
    dUmax_(ptf.dUmax_),
    Umin_(ptf.Umin_),
    Umax_(ptf.Umax_),
    updateEveryNTimeSteps_(ptf.updateEveryNTimeSteps_),
    lastUpdateTimeIndex_(ptf.lastUpdateTimeIndex_),
    initFromPatchAverage_(ptf.initFromPatchAverage_),
    initialized_(ptf.initialized_),
    U_(ptf.U_),
    eInt_(ptf.eInt_),
    verbose_(ptf.verbose_),
    lastPrintTimeIndex_(ptf.lastPrintTimeIndex_),
    useExternalDict_(ptf.useExternalDict_),
    externalDictName_(ptf.externalDictName_),
    extDictPtr_(nullptr),
    extDictEventNo_(-1)          // IMPORTANT
{}

void Foam::pidUNeFvPatchScalarField::updateCoeffs()
{
    if (updated()) return;

    const fvMesh& m = patch().boundaryMesh().mesh();
    const label ti = m.time().timeIndex();

    // 1) Always read/apply external dict (does NOT update U_ by itself)
    //    (applyDict must NOT touch lastUpdateTimeIndex_ unless resetController==true)
    readAndApplyExternalDict();

    // 2) Compute cadence for controller update (this controls when U_ changes)
    bool doUpdate = false;
    if (lastUpdateTimeIndex_ < 0 || (ti - lastUpdateTimeIndex_) >= updateEveryNTimeSteps_)
    {
        doUpdate = true;
        // IMPORTANT: only set this when we actually update U_
        lastUpdateTimeIndex_ = ti;
    }

    // 3) One-time initialization from patchInternalField (optional)
    //    (This is independent of doUpdate; it sets initial U_ once)
    if (initFromPatchAverage_ && !initialized_)
    {
        tmp<scalarField> tUi = this->patchInternalField();
        const scalarField& Ui = tUi();

        scalar sumU = 0.0, nU = 0.0;
        forAll(Ui, i) { sumU += Ui[i]; nU += 1.0; }

        reduce(sumU, sumOp<scalar>());
        reduce(nU,  sumOp<scalar>());

        if (nU > SMALL) U_ = sumU/nU;

        U_ = max(Umin_, min(Umax_, U_));
        initialized_ = true;

        if (verbose_ && Pstream::master())
        {
            Info<< "pidUNe: init U from patchInternal avg = " << U_ << nl;
        }
    }

    // 4) Only update controller (and thus U_) on cadence
    scalar Iabs = 0.0, A = 0.0, Igoal = 0.0, jAvg = 0.0, e = 0.0;

    if (doUpdate)
    {
        // Measure current through patch
        if (!m.foundObject<volVectorField>(currentFieldName_))
        {
            FatalErrorInFunction
                << "Field '" << currentFieldName_ << "' not found as volVectorField."
                << exit(FatalError);
        }

        const volVectorField& J = m.lookupObject<volVectorField>(currentFieldName_);
        const surfaceVectorField Jf(fvc::interpolate(J));

        const label patchI = patch().index();
        const vectorField& Sf  = patch().Sf();
        const vectorField& Jfp = Jf.boundaryField()[patchI];

        scalar I_local = 0.0, A_local = 0.0;
        forAll(Sf, i)
        {
            I_local += (Jfp[i] & Sf[i]);
            A_local += mag(Sf[i]);
        }

        reduce(I_local, sumOp<scalar>());
        reduce(A_local, sumOp<scalar>());

        Iabs = mag(I_local);
        A    = A_local;

        // Error signal
        if (mode_ == cmTotal)
        {
            Igoal = jGoal_*A;
            e = Igoal - Iabs;
        }
        else // cmAverage
        {
            if (A > SMALL) jAvg = Iabs/A;
            e = jGoal_ - jAvg;
        }

        // PI update (only here!)
        if (A > SMALL)
        {
            const scalar dt = m.time().deltaTValue();

            // anti-windup
            if (!((U_ <= Umin_ && e < 0) || (U_ >= Umax_ && e > 0)))
            {
                eInt_ += e*dt;
            }

            scalar dU = Kp_*e + Ki_*eInt_;
            dU = max(-dUmax_, min(dUmax_, dU));

            U_ -= dU;
            U_ = max(Umin_, min(Umax_, U_));
        }

        // Optional: print only on updates (recommended)
        if (Pstream::master())
        {
            if (mode_ == cmTotal)
            {
                Info<< "pidUNe["<< controlModeNames_[mode_] << "] patch="<< patch().name()
                    << " t="<< m.time().timeName()
                    << " I="<< Iabs << " A  A="<< A << " m2  Igoal="<< Igoal
                    << " e="<< e << " U="<< U_ << nl;
            }
            else
            {
                Info<< "pidUNe["<< controlModeNames_[mode_] << "] patch="<< patch().name()
                    << " t="<< m.time().timeName()
                    << " jAvg="<< jAvg << " A/m2  A="<< A << " m2  jGoal="<< jGoal_
                    << " e="<< e << " U="<< U_ << nl;
            }
                // New line: PID settings
        Info<< "  PID: Kp=" << Kp_
        << " Ki=" << Ki_
        << " dUmax=" << dUmax_
        << " Umin=" << Umin_
        << " Umax=" << Umax_
        << " updateEveryNTimeSteps=" << updateEveryNTimeSteps_
        << " currentField=" << currentFieldName_
        << " eInt=" << eInt_
        << nl;
        }
    }

    // 5) Always set patch to current U_ (even when not updating)
    this->operator==(U_);

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::pidUNeFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);

    os.writeEntry("controlMode", controlModeNames_[mode_]);
    os.writeEntry("currentField", currentFieldName_);
    os.writeEntry("jGoal", jGoal_);
    os.writeEntry("Kp", Kp_);
    os.writeEntry("Ki", Ki_);
    os.writeEntry("dUmax", dUmax_);
    os.writeEntry("Umin", Umin_);
    os.writeEntry("Umax", Umax_);
    os.writeEntry("initialU", U_);
    os.writeEntry("updateEveryNTimeSteps", updateEveryNTimeSteps_);
    os.writeEntry("initFromPatchAverage", initFromPatchAverage_);
    os.writeEntry("verbose", verbose_);

    // external dict settings (so it’s reproducible)
    os.writeEntry("useExternalDict", useExternalDict_);
    os.writeEntry("externalDictName", externalDictName_);

    writeEntry("value", os);
}

// ************************************************************************* //
