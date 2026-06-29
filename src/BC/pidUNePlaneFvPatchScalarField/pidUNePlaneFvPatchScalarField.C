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

\*---------------------------------------------------------------------------*/

#include "pidUNePlaneFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvCFD.H"
#include "faceZoneMesh.H"
#include "faceZone.H"
#include "fvc.H"
#include "PstreamReduceOps.H"

namespace Foam
{
    makePatchTypeField(fvPatchScalarField, pidUNePlaneFvPatchScalarField);
}

const Foam::Enum<Foam::pidUNePlaneFvPatchScalarField::controlMode>
Foam::pidUNePlaneFvPatchScalarField::controlModeNames_
({
    { controlMode::cmTotal,   "total"   },
    { controlMode::cmAverage, "average" }
});


void Foam::pidUNePlaneFvPatchScalarField::applyDict(const dictionary& d)
{
    mode_ = controlModeNames_.getOrDefault("controlMode", d, mode_);

    d.readIfPresent("currentField", currentFieldName_);
    d.readIfPresent("faceZoneName", faceZoneName_);
    d.readIfPresent("useProjected", useProjected_);

    if (d.readIfPresent("normal", nHat_))
    {
        if (mag(nHat_) > VSMALL)
        {
            nHat_ /= mag(nHat_);
        }
        else
        {
            nHat_ = vector(0, 0, 1);
        }
    }

    d.readIfPresent("jGoal", jGoal_);
    d.readIfPresent("Kp", Kp_);
    d.readIfPresent("Ki", Ki_);
    d.readIfPresent("dUmax", dUmax_);
    d.readIfPresent("Umin", Umin_);
    d.readIfPresent("Umax", Umax_);
    d.readIfPresent("updateEveryNTimeSteps", updateEveryNTimeSteps_);
    d.readIfPresent("initFromPatchAverage", initFromPatchAverage_);
    d.readIfPresent("verbose", verbose_);

    if (updateEveryNTimeSteps_ < 1)
    {
        updateEveryNTimeSteps_ = 1;
    }

    bool resetController(false);
    if (d.readIfPresent("resetController", resetController) && resetController)
    {
        d.readIfPresent("initialU", U_);
        eInt_ = 0.0;
        initialized_ = false;
        firstControllerUpdateSkipped_ = false;
        lastUpdateTimeIndex_ = -1;
    }

    U_ = max(Umin_, min(Umax_, U_));
}


const Foam::dictionary*
Foam::pidUNePlaneFvPatchScalarField::externalDictTop() const
{
    if (!useExternalDict_)
    {
        return nullptr;
    }

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


void Foam::pidUNePlaneFvPatchScalarField::ensureExternalDict()
{
    if (!useExternalDict_)
    {
        return;
    }

    const fvMesh& m = patch().boundaryMesh().mesh();

    if (m.foundObject<IOdictionary>(externalDictName_))
    {
        return;
    }

    if (extDictPtr_.valid())
    {
        return;
    }

    extDictPtr_.reset
    (
        new IOdictionary
        (
            IOobject
            (
                externalDictName_,
                m.time().system(),
                m,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE
            )
        )
    );
}


void Foam::pidUNePlaneFvPatchScalarField::readAndApplyExternalDict()
{
    if (!useExternalDict_)
    {
        return;
    }

    ensureExternalDict();

    const fvMesh& m = patch().boundaryMesh().mesh();
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
    if (!topPtr)
    {
        return;
    }

    const dictionary& top = *topPtr;
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

    applyDict(*usePtr);

    if (verbose_ && Pstream::master())
    {
        const scalar jFile = usePtr->lookupOrDefault<scalar>("jGoal", jGoal_);

        Info<< "pidUNePlane: applied " << externalDictName_
            << " patch=" << patch().name()
            << " jGoal(file)=" << jFile
            << " jGoal(mem)=" << jGoal_
            << " faceZone=" << faceZoneName_
            << " useProjected=" << useProjected_
            << " normal=" << nHat_
            << " updateN=" << updateEveryNTimeSteps_
            << " path=\"" << dictPath << "\""
            << nl;
    }
}


void Foam::pidUNePlaneFvPatchScalarField::calculatePlaneCurrent
(
    scalar& Iabs,
    scalar& A,
    scalar& jAvg
) const
{
    const fvMesh& m = patch().boundaryMesh().mesh();

    Iabs = 0.0;
    A = 0.0;
    jAvg = 0.0;

    if (!m.foundObject<volVectorField>(currentFieldName_))
    {
        FatalErrorInFunction
            << "Field '" << currentFieldName_
            << "' not found as volVectorField."
            << exit(FatalError);
    }

    const volVectorField& J = m.lookupObject<volVectorField>(currentFieldName_);

    const faceZoneMesh& fzMesh = m.faceZones();
    const label fzId = fzMesh.findZoneID(faceZoneName_);

    if (fzId < 0)
    {
        FatalErrorInFunction
            << "faceZone '" << faceZoneName_ << "' not found."
            << exit(FatalError);
    }

    const faceZone& fz = fzMesh[fzId];

    const labelList faces(fz);
    const boolList  flip(fz.flipMap());

    if (flip.size() != faces.size())
    {
        FatalErrorInFunction
            << "faceZone '" << faceZoneName_
            << "': flipMap size " << flip.size()
            << " != faces size " << faces.size()
            << exit(FatalError);
    }

    const surfaceVectorField Jf(fvc::interpolate(J));
    const surfaceVectorField& Sf = m.Sf();

    if
    (
        Sf.primitiveField().size() != m.nInternalFaces()
     || Jf.primitiveField().size() != m.nInternalFaces()
    )
    {
        FatalErrorInFunction
            << "Unexpected primitive sizes: "
            << "Sf.primitiveField=" << Sf.primitiveField().size()
            << " Jf.primitiveField=" << Jf.primitiveField().size()
            << " nInternalFaces=" << m.nInternalFaces()
            << exit(FatalError);
    }

    auto SfAtFace = [&](const label faceI) -> vector
    {
        if (faceI < m.nInternalFaces())
        {
            return Sf.primitiveField()[faceI];
        }

        const label patchI = m.boundaryMesh().whichPatch(faceI);
        if (patchI < 0)
        {
            FatalErrorInFunction
                << "Boundary face " << faceI
                << " has no patch (whichPatch=-1)"
                << exit(FatalError);
        }

        const polyPatch& pp = m.boundaryMesh()[patchI];
        const label local = faceI - pp.start();

        if (local < 0 || local >= pp.size())
        {
            FatalErrorInFunction
                << "Bad boundary local index: face=" << faceI
                << " patch=" << pp.name()
                << " start=" << pp.start()
                << " size=" << pp.size()
                << " local=" << local
                << exit(FatalError);
        }

        return Sf.boundaryField()[patchI][local];
    };

    auto JfAtFace = [&](const label faceI) -> vector
    {
        if (faceI < m.nInternalFaces())
        {
            return Jf.primitiveField()[faceI];
        }

        const label patchI = m.boundaryMesh().whichPatch(faceI);
        if (patchI < 0)
        {
            FatalErrorInFunction
                << "Boundary face " << faceI
                << " has no patch (whichPatch=-1)"
                << exit(FatalError);
        }

        const polyPatch& pp = m.boundaryMesh()[patchI];
        const label local = faceI - pp.start();

        if (local < 0 || local >= pp.size())
        {
            FatalErrorInFunction
                << "Bad boundary local index: face=" << faceI
                << " patch=" << pp.name()
                << " start=" << pp.start()
                << " size=" << pp.size()
                << " local=" << local
                << exit(FatalError);
        }

        return Jf.boundaryField()[patchI][local];
    };

    scalar I_local = 0.0;
    scalar A_local = 0.0;

    forAll(faces, i)
    {
        const label faceI = faces[i];

        vector s = SfAtFace(faceI);
        if (flip[i])
        {
            s = -s;
        }

        const vector jf = JfAtFace(faceI);

        if (useProjected_)
        {
            // Plane flux: (J.nHat)*dAplane with dAplane = Sf.nHat.
            const scalar dA = (s & nHat_);
            I_local += (jf & nHat_)*dA;
            A_local += mag(dA);
        }
        else
        {
            // Actual mesh-face flux through the faceZone faces.
            I_local += (jf & s);
            A_local += mag(s);
        }
    }

    reduce(I_local, sumOp<scalar>());
    reduce(A_local, sumOp<scalar>());

    Iabs = mag(I_local);
    A = A_local;

    if (A > SMALL)
    {
        jAvg = Iabs/A;
    }
}


Foam::pidUNePlaneFvPatchScalarField::pidUNePlaneFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    mode_(cmTotal),
    currentFieldName_("Ie"),
    faceZoneName_("midPlaneZone"),
    useProjected_(true),
    nHat_(vector(0, 0, 1)),
    jGoal_(0.0),
    Kp_(0.0),
    Ki_(0.0),
    dUmax_(0.1),
    Umin_(-GREAT),
    Umax_(GREAT),
    updateEveryNTimeSteps_(1),
    firstControllerUpdateSkipped_(false),
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
    extDictEventNo_(-1)
{}


Foam::pidUNePlaneFvPatchScalarField::pidUNePlaneFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    mode_(controlModeNames_.getOrDefault("controlMode", dict, cmTotal)),
    currentFieldName_(dict.getOrDefault<word>("currentField", "Ie")),
    faceZoneName_(dict.getOrDefault<word>("faceZoneName", "midPlaneZone")),
    useProjected_(dict.getOrDefault<bool>("useProjected", true)),
    nHat_(dict.getOrDefault<vector>("normal", vector(0, 0, 1))),
    jGoal_(dict.getOrDefault<scalar>("jGoal", 0.0)),
    Kp_(dict.getOrDefault<scalar>("Kp", 0.0)),
    Ki_(dict.getOrDefault<scalar>("Ki", 0.0)),
    dUmax_(dict.getOrDefault<scalar>("dUmax", 0.1)),
    Umin_(dict.getOrDefault<scalar>("Umin", -GREAT)),
    Umax_(dict.getOrDefault<scalar>("Umax", GREAT)),
    updateEveryNTimeSteps_(dict.getOrDefault<label>("updateEveryNTimeSteps", 1)),
    firstControllerUpdateSkipped_
    (
        dict.getOrDefault<bool>("firstControllerUpdateSkip", true)
    ),
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
    extDictEventNo_(-1)
{
    if (mag(nHat_) > VSMALL)
    {
        nHat_ /= mag(nHat_);
    }
    else
    {
        nHat_ = vector(0, 0, 1);
    }

    if (updateEveryNTimeSteps_ < 1)
    {
        updateEveryNTimeSteps_ = 1;
    }

    if (!dict.found("value"))
    {
        this->operator==(U_);
    }
}


Foam::pidUNePlaneFvPatchScalarField::pidUNePlaneFvPatchScalarField
(
    const pidUNePlaneFvPatchScalarField& ptf
)
:
    fixedValueFvPatchScalarField(ptf),
    mode_(ptf.mode_),
    currentFieldName_(ptf.currentFieldName_),
    faceZoneName_(ptf.faceZoneName_),
    useProjected_(ptf.useProjected_),
    nHat_(ptf.nHat_),
    jGoal_(ptf.jGoal_),
    Kp_(ptf.Kp_),
    Ki_(ptf.Ki_),
    dUmax_(ptf.dUmax_),
    Umin_(ptf.Umin_),
    Umax_(ptf.Umax_),
    updateEveryNTimeSteps_(ptf.updateEveryNTimeSteps_),
    firstControllerUpdateSkipped_(ptf.firstControllerUpdateSkipped_),
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
    extDictEventNo_(-1)
{}


Foam::pidUNePlaneFvPatchScalarField::pidUNePlaneFvPatchScalarField
(
    const pidUNePlaneFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    mode_(ptf.mode_),
    currentFieldName_(ptf.currentFieldName_),
    faceZoneName_(ptf.faceZoneName_),
    useProjected_(ptf.useProjected_),
    nHat_(ptf.nHat_),
    jGoal_(ptf.jGoal_),
    Kp_(ptf.Kp_),
    Ki_(ptf.Ki_),
    dUmax_(ptf.dUmax_),
    Umin_(ptf.Umin_),
    Umax_(ptf.Umax_),
    updateEveryNTimeSteps_(ptf.updateEveryNTimeSteps_),
    firstControllerUpdateSkipped_(ptf.firstControllerUpdateSkipped_),
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
    extDictEventNo_(-1)
{}


Foam::pidUNePlaneFvPatchScalarField::pidUNePlaneFvPatchScalarField
(
    const pidUNePlaneFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    mode_(ptf.mode_),
    currentFieldName_(ptf.currentFieldName_),
    faceZoneName_(ptf.faceZoneName_),
    useProjected_(ptf.useProjected_),
    nHat_(ptf.nHat_),
    jGoal_(ptf.jGoal_),
    Kp_(ptf.Kp_),
    Ki_(ptf.Ki_),
    dUmax_(ptf.dUmax_),
    Umin_(ptf.Umin_),
    Umax_(ptf.Umax_),
    updateEveryNTimeSteps_(ptf.updateEveryNTimeSteps_),
    firstControllerUpdateSkipped_(ptf.firstControllerUpdateSkipped_),
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
    extDictEventNo_(-1)
{}


void Foam::pidUNePlaneFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvMesh& m = patch().boundaryMesh().mesh();
    const label ti = m.time().timeIndex();

    // Always read/apply external dict.  This does not update U_ unless
    // resetController is true.
    readAndApplyExternalDict();

    bool doUpdate = false;

    if (lastUpdateTimeIndex_ < 0)
    {
        // First call: do not run PID because currentField may not be valid yet.
        lastUpdateTimeIndex_ = ti;
        firstControllerUpdateSkipped_ = true;

        if (verbose_ && Pstream::master())
        {
            Info<< "pidUNePlane: skipping first PID update at timeIndex="
                << ti << " because " << currentFieldName_
                << " may not be calculated yet" << nl;
        }
    }
    else if ((ti - lastUpdateTimeIndex_) >= updateEveryNTimeSteps_)
    {
        doUpdate = true;
        lastUpdateTimeIndex_ = ti;
    }

    if (initFromPatchAverage_ && !initialized_)
    {
        tmp<scalarField> tUi = this->patchInternalField();
        const scalarField& Ui = tUi();

        scalar sumU = 0.0;
        scalar nU = 0.0;

        forAll(Ui, i)
        {
            sumU += Ui[i];
            nU += 1.0;
        }

        reduce(sumU, sumOp<scalar>());
        reduce(nU,  sumOp<scalar>());

        if (nU > SMALL)
        {
            U_ = sumU/nU;
        }

        U_ = max(Umin_, min(Umax_, U_));
        initialized_ = true;

        if (verbose_ && Pstream::master())
        {
            Info<< "pidUNePlane: init U from patchInternal avg = "
                << U_ << nl;
        }
    }

    scalar Iabs = 0.0;
    scalar A = 0.0;
    scalar Igoal = 0.0;
    scalar jAvg = 0.0;
    scalar e = 0.0;

    if (doUpdate)
    {
        calculatePlaneCurrent(Iabs, A, jAvg);

        if (mode_ == cmTotal)
        {
            Igoal = jGoal_*A;
            e = Igoal - Iabs;
        }
        else
        {
            e = jGoal_ - jAvg;
        }

        if (A > SMALL)
        {
            const scalar dt = m.time().deltaTValue();

            // Anti-windup.
            if (!((U_ <= Umin_ && e < 0) || (U_ >= Umax_ && e > 0)))
            {
                eInt_ += e*dt;
            }

            scalar dU = Kp_*e + Ki_*eInt_;
            dU = max(-dUmax_, min(dUmax_, dU));

            U_ -= dU;
            U_ = max(Umin_, min(Umax_, U_));
        }

        if (verbose_ && Pstream::master())
        {
            if (mode_ == cmTotal)
            {
                Info<< "pidUNePlane[" << controlModeNames_[mode_]
                    << "] patch=" << patch().name()
                    << " faceZone=" << faceZoneName_
                    << " t=" << m.time().timeName()
                    << " I=" << Iabs << " A  A=" << A
                    << " m2  Igoal=" << Igoal
                    << " e=" << e << " U=" << U_ << nl;
            }
            else
            {
                Info<< "pidUNePlane[" << controlModeNames_[mode_]
                    << "] patch=" << patch().name()
                    << " faceZone=" << faceZoneName_
                    << " t=" << m.time().timeName()
                    << " jAvg=" << jAvg << " A/m2  A=" << A
                    << " m2  jGoal=" << jGoal_
                    << " e=" << e << " U=" << U_ << nl;
            }

            Info<< "  PID: Kp=" << Kp_
                << " Ki=" << Ki_
                << " dUmax=" << dUmax_
                << " Umin=" << Umin_
                << " Umax=" << Umax_
                << " updateEveryNTimeSteps=" << updateEveryNTimeSteps_
                << " currentField=" << currentFieldName_
                << " useProjected=" << useProjected_
                << " normal=" << nHat_
                << " eInt=" << eInt_
                << nl;
        }
    }

    this->operator==(U_);

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::pidUNePlaneFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);

    os.writeEntry("controlMode", controlModeNames_[mode_]);
    os.writeEntry("currentField", currentFieldName_);
    os.writeEntry("faceZoneName", faceZoneName_);
    os.writeEntry("useProjected", useProjected_);
    os.writeEntry("normal", nHat_);
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
    os.writeEntry("useExternalDict", useExternalDict_);
    os.writeEntry("externalDictName", externalDictName_);

    writeEntry("value", os);
}

// ************************************************************************* //
