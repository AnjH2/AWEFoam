/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
\*---------------------------------------------------------------------------*/

// ================= pidUNeConstraint.C =================
#include "pidUNeConstraint.H"
#include "addToRunTimeSelectionTable.H"
#include "fvCFD.H"
#include "faceZoneMesh.H"
#include "faceZone.H"
#include "fvc.H"
#include "PstreamReduceOps.H"

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(pidUNeConstraint, 0);
    addToRunTimeSelectionTable(option, pidUNeConstraint, dictionary);
}
}

// --- Enum mapping must be at file scope (not inside a function)
const Foam::Enum<Foam::fv::pidUNeConstraint::controlMode>
Foam::fv::pidUNeConstraint::controlModeNames_
({
    { controlMode::cmTotalCurrent, "total" },
    { controlMode::cmAverageFlux,  "average" }
});


Foam::fv::pidUNeConstraint::pidUNeConstraint
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fv::cellSetOption(name, modelType, dict, mesh),

    mode_(controlModeNames_.getOrDefault("controlMode", coeffs_, cmTotalCurrent)),
    useProjected_(coeffs_.getOrDefault<bool>("useProjected", true)),
    nHat_(vector(0,0,1)),

    faceZoneName_(coeffs_.getOrDefault<word>("faceZoneName", "midPlaneZone")),
    currentFieldName_(coeffs_.getOrDefault<word>("currentField", "Ie")),

    jGoal_(coeffs_.getOrDefault<scalar>("jGoal", 3000)),
    Kp_(coeffs_.getOrDefault<scalar>("Kp", 1e-2)),
    Ki_(coeffs_.getOrDefault<scalar>("Ki", 0.0)),
    dUmax_(coeffs_.getOrDefault<scalar>("dUmax", 0.1)),
    Umin_(coeffs_.getOrDefault<scalar>("Umin", -10.0)),
    Umax_(coeffs_.getOrDefault<scalar>("Umax", -1.35)),
    initialU_(coeffs_.getOrDefault<scalar>("initialU", -1.35)),
    Nupdate_(coeffs_.getOrDefault<label>("Nupdate", 8)),

    verbose_(coeffs_.getOrDefault<bool>("verbose", true)),
    printEveryTimeStep_(coeffs_.getOrDefault<bool>("printEveryTimeStep", true)),
    printOnUpdate_(coeffs_.getOrDefault<bool>("printOnUpdate", true)),
    lastPrintTimeIndex_(-1),

    U_(initialU_),
    eInt_(0.0),
    k_(0),
    
    initFromField_(coeffs_.getOrDefault<bool>("initFromField", false)),
    initialized_(false),
    
    updateEveryNTimeSteps_(coeffs_.getOrDefault<label>("updateEveryNTimeSteps", 1)),
    lastUpdateTimeIndex_(-1)
    
    
    
{
    if (Nupdate_ < 1) Nupdate_ = 1;

    // allow overriding normal in coeffs
    coeffs_.readIfPresent("normal", nHat_);
    if (mag(nHat_) > VSMALL) nHat_ /= mag(nHat_);
    else nHat_ = vector(0,0,1);

    // target fields (needed so constrain() is called)
    wordList flds;
    if (dict.readIfPresent("fields", flds))
    {
        fieldNames_ = flds;
    }
    else if (coeffs_.readIfPresent("fields", flds))
    {
        fieldNames_ = flds;
    }
    else
    {
        fieldNames_.resize(1);
        fieldNames_[0] = "UNe";
    }

    if (Pstream::master())
    {
        Info<< "pidUNeConstraint constructed: name=" << name
            << " type=" << modelType << nl
            << "  targetFields=" << fieldNames_ << nl
            << "  faceZone=" << faceZoneName_
            << " currentField=" << currentFieldName_ << nl
            << "  nCellsSelected=" << cells_.size() << nl;
    }

    fv::option::resetApplied();
    
    if (updateEveryNTimeSteps_ < 1) updateEveryNTimeSteps_ = 1;
}


bool Foam::fv::pidUNeConstraint::read(const dictionary& dict)
{
    if (!fv::cellSetOption::read(dict)) return false;

    coeffs_.readIfPresent("faceZoneName", faceZoneName_);
    coeffs_.readIfPresent("currentField", currentFieldName_);

    coeffs_.readIfPresent("jGoal", jGoal_);
    coeffs_.readIfPresent("Kp", Kp_);
    coeffs_.readIfPresent("Ki", Ki_);
    coeffs_.readIfPresent("dUmax", dUmax_);
    coeffs_.readIfPresent("Umin", Umin_);
    coeffs_.readIfPresent("Umax", Umax_);
    coeffs_.readIfPresent("initialU", initialU_);
    coeffs_.readIfPresent("Nupdate", Nupdate_);

    coeffs_.readIfPresent("verbose", verbose_);
    coeffs_.readIfPresent("printEveryTimeStep", printEveryTimeStep_);
    coeffs_.readIfPresent("printOnUpdate", printOnUpdate_);

    mode_ = controlModeNames_.getOrDefault("controlMode", coeffs_, cmTotalCurrent);
    coeffs_.readIfPresent("useProjected", useProjected_);
    
    coeffs_.readIfPresent("updateEveryNTimeSteps", updateEveryNTimeSteps_);
    if (updateEveryNTimeSteps_ < 1) updateEveryNTimeSteps_ = 1;
    

    if (coeffs_.readIfPresent("normal", nHat_))
    {
        if (mag(nHat_) > VSMALL) nHat_ /= mag(nHat_);
        else nHat_ = vector(0,0,1);
    }

    if (Nupdate_ < 1) Nupdate_ = 1;

    wordList flds;
    if (dict.readIfPresent("fields", flds))
    {
        fieldNames_ = flds;
    }
    else if (coeffs_.readIfPresent("fields", flds))
    {
        fieldNames_ = flds;
    }
    else
    {
        fieldNames_.resize(1);
        fieldNames_[0] = "UNe";
    }

    fv::option::resetApplied();
    return true;
}


void Foam::fv::pidUNeConstraint::constrain(fvMatrix<scalar>& eqn, const label)
{
    const fvMesh& m = mesh();
    const label ti = m.time().timeIndex();
    if (initFromField_ && !initialized_)
    {
        if (!m.foundObject<volScalarField>("UNe"))
        {
         FatalErrorInFunction
             << "initFieldName '" << "UNe" << "' not found as volScalarField."
             << exit(FatalError);
        }

        const volScalarField& Ufld = m.lookupObject<volScalarField>("UNe");

        // Average over selected cells
        scalar sumU = 0.0;
        scalar nU = 0.0;
        forAll(cells_, i)
        {
            const label c = cells_[i];
            sumU += Ufld[c];
            nU += 1.0;
        }
        reduce(sumU, sumOp<scalar>());
        reduce(nU, sumOp<scalar>());

        if (nU > SMALL)
        {
            U_ = sumU/nU;
            U_ = max(Umin_, min(Umax_, U_));
        }

        initialized_ = true;

        if (verbose_ && Pstream::master())
        {
            Info<< "pidUNeConstraint: initialized UNe from field '"
                << "UNe" << "' avg over cellZone = " << U_ << nl;
        }
    }
    // --- current density field (cell)
    if (!m.foundObject<volVectorField>(currentFieldName_))
    {
        FatalErrorInFunction
            << "Field '" << currentFieldName_
            << "' not found as volVectorField."
            << exit(FatalError);
    }
    const volVectorField& J = m.lookupObject<volVectorField>(currentFieldName_);

    // --- faceZone
    const faceZoneMesh& fzMesh = m.faceZones();
    const label fzId = fzMesh.findZoneID(faceZoneName_);
    if (fzId < 0)
    {
        FatalErrorInFunction
            << "faceZone '" << faceZoneName_ << "' not found."
            << exit(FatalError);
    }
    const faceZone& fz = fzMesh[fzId];

    // Copy lists (avoid any lifetime issues)
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

    // Interpolate J to faces (surface field with internal+boundary storage)
    const surfaceVectorField Jf(fvc::interpolate(J));
    const surfaceVectorField& Sf = m.Sf();

    // Sanity: primitive (internal) sizes should match nInternalFaces
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

    // Helper: get surface field value by GLOBAL face label
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
                << "Boundary face " << faceI << " has no patch (whichPatch=-1)"
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
                << "Boundary face " << faceI << " has no patch (whichPatch=-1)"
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

    // --- integrate current and area over faceZone

scalar I_local = 0.0;
scalar A_local = 0.0;

forAll(faces, i)
{
    const label faceI = faces[i];

    vector s = SfAtFace(faceI);
    if (flip[i]) s = -s;

    const vector jf = JfAtFace(faceI);

    if (useProjected_)
    {
        // plane flux: (J·nHat) dA_plane where dA_plane = (Sf·nHat)
        const scalar dA = (s & nHat_);     // signed projected area
        I_local += (jf & nHat_) * dA;      // [A]
        A_local += mag(dA);                // [m2] projected
    }
    else
    {
        // face-normal flux through the actual mesh faces
        I_local += (jf & s);               // [A]
        A_local += mag(s);                 // [m2]
    }
}

reduce(I_local, sumOp<scalar>());
reduce(A_local, sumOp<scalar>());

// Measured values
const scalar I = mag(I_local);     // keep sign; use mag(I) if you want magnitude
const scalar A = A_local;

// Compute error depending on mode
scalar e = 0.0;
scalar Igoal = 0.0;
scalar jAvg = 0.0;

if (mode_ == cmTotalCurrent)
{
    Igoal = jGoal_*A;         // [A]
    e = Igoal - I;            // [A]
}
else // cmAverageFlux
{
    if (A > SMALL) jAvg = I/A; // [A/m2]
    e = jGoal_ - jAvg;         // [A/m2]
}

// PI update uses e (units depend on mode)
const scalar dt = m.time().deltaTValue();

bool doUpdate = false;
if (A > SMALL)
{
    // Update if first time, or enough timesteps since last update
    if (lastUpdateTimeIndex_ < 0 || (ti - lastUpdateTimeIndex_) >= updateEveryNTimeSteps_)
    {
        doUpdate = true;
        lastUpdateTimeIndex_ = ti;
    }
}

if (doUpdate)
{
    if (!((U_ <= Umin_ && e < 0) || (U_ >= Umax_ && e > 0)))
    {
        eInt_ += e*dt;
    }

    scalar dU = Kp_*e + Ki_*eInt_;
    dU = max(-dUmax_, min(dUmax_, dU));

    U_ -= dU;
    U_ = max(Umin_, min(Umax_, U_));
}
else
{
    U_ = max(Umin_, min(Umax_, U_));
}

// Print (master only)
if (verbose_ && Pstream::master() && m.time().timeIndex() != lastPrintTimeIndex_)
{
    lastPrintTimeIndex_ = m.time().timeIndex();

    if (mode_ == cmTotalCurrent)
    {
        Info<< "pidUNeConstraint[" << controlModeNames_[mode_] << "] (t="<< m.time().timeName() << "): "
            << "I="<< I << " A, A="<< A << " m2, Igoal="<< Igoal << " A, e="<< e
            << ", eInt="<< eInt_ << ", U="<< U_ << " V\n";
    }
    else
    {
        Info<< "pidUNeConstraint[" << controlModeNames_[mode_] << "] (t="<< m.time().timeName() << "): "
            << "jAvg="<< jAvg << " A/m2, A="<< A << " m2, jGoal="<< jGoal_ << " A/m2, e="<< e
            << ", eInt="<< eInt_ << ", U="<< U_ << " V\n";
    }
}


    // --- apply constraint
    eqn.setValues(cells_, U_);
}

