Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading porousProperties\n" << endl;
porousProperties poroM(mesh);

#include "createPhi.H"


Info<< "Reading incompressibleTwoPhaseInteractingMixture\n" << endl;
incompressibleTwoPhaseInteractingMixture mixture(U, phi);

Info<< "Reading reactionProperties\n" << endl;
reactionProperties react(mesh);

volScalarField& alpha1(mixture.alpha1());

// Mixture density
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mixture.rho()
);

// Mass flux
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);

surfaceScalarField alpha2Phi
(
    IOobject
    (
        "alpha2Phi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    phi
);

surfaceScalarField alpha1Phi
(
    IOobject
    (
        "alpha1Phi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    phi
);

// Relative Velocity
autoPtr<relativeVelocityModel> UdmModelPtr
(
    relativeVelocityModel::New
    (
        mixture,
        mixture
    )
);


// Construct compressible turbulence model
autoPtr
<
    CompressibleTurbulenceModel<incompressibleTwoPhaseInteractingMixture>
> turbulence
(
    CompressibleTurbulenceModel<incompressibleTwoPhaseInteractingMixture>
    ::New(rho, U, rhoPhi, mixture)
);


#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"


volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
    p_rgh = p - rho*gh;
}

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());

// MULES Correction
tmp<surfaceScalarField> talphaPhiCorr0;

//********** Potential related Fields **********//

Info << "Reading potential field UPe \n" << endl;
volScalarField UPe
(
	IOobject
	(
		"UPe",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
	),
	mesh
);

Info << "Reading potential field UNe \n" << endl;
volScalarField UNe
(
	IOobject
	(
		"UNe",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
	),
	mesh
);
Info << "Reading potential field Ue \n" << endl;
volScalarField Ue
(
	IOobject
	(
		"Ue",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
	),
	mesh
);
Info << "Reading or preparing transfer current field J \n" << endl;
volScalarField J
(
	IOobject
	(
		"J",
		runTime.timeName(),
		mesh,
		IOobject::READ_IF_PRESENT,
		IOobject::AUTO_WRITE
	),
	mesh,
	dimensionSet(0,-3,0,0,0,1,0)
);

Info << "Reading or preparing overpotential field eta \n" << endl;
volScalarField eta
(
	IOobject
	(
		"eta",
		runTime.timeName(),
		mesh,
		IOobject::READ_IF_PRESENT,
		IOobject::AUTO_WRITE
	),
	mesh,
	dimensionSet(1,2,-3,0,0,-1,0)
);

#include "createMRF.H"
#include "createFvOptions.H"
