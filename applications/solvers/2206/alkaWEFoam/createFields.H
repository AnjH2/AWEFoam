// Create the primary hydrodynamic, electrochemical and transport fields and
// construct the runtime-selectable submodels used by alkaWEFoam.
IOdictionary transportProperties
(
	IOobject
	(
		"transportProperties",
		runTime.constant(),
		mesh,
		IOobject::MUST_READ_IF_MODIFIED,
		IOobject::NO_WRITE
	)
);
IOdictionary fvSolution
(
	IOobject
	(
		"fvSolution",
		runTime.system(),
		mesh,
		IOobject::MUST_READ_IF_MODIFIED,
		IOobject::NO_WRITE
	)
);



//scalar to define number of iterations each variables will be solved

const scalar speciesCorrections=fvSolution.subDict("alkaWEFoam").get<scalar>("speciesCorrections");
const scalar potentialCorrections=fvSolution.subDict("alkaWEFoam").get<scalar>("potentialCorrections");
const scalar outerChemicalCorrections=fvSolution.subDict("alkaWEFoam").get<scalar>("outerChemicalCorrections");
//scalar nT=transportProperties.get<scalar>("n_Temp");
Info<< "Reading field p_rgh\n" << endl;

volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
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
    p_rgh
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
volVectorField Ur
(
    IOobject
    (
        "Ur",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    U
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
volScalarField Solid
(
    IOobject
    (
        "Solid",
        mesh.time().timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar ("Solid", dimless, 1)
);
Info<< "Reading porousProperties\n" << endl;
porousProperties poroM(mesh);

Info<< "Reading speciesProperties\n" << endl;
speciesProperties speciesP(mesh,transportProperties);
#include "createPhi.H"

volScalarField divU
(
    IOobject
    (
        "divU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fvc::div(phi)
);
Info<< "Reading incompressibleTwoPhaseInteractingMixture\n" << endl;
incompressibleTwoPhaseInteractingMixture mixture(U, phi,transportProperties);

Info<< "Reading reactionProperties\n" << endl;
reactionProperties react(mesh);



volScalarField& alpha1(mixture.alpha1());

// Mixture density
Info<<"Creating field density"<<endl;
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mixture.rho()*poroM.eps()
);

// Mass flux
Info<<"Creating mass flux rhoPhi"<<endl;
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
    fvc::interpolate(mixture.rho())*phi
);

 
Info<<"Creating Saturation flux"<<endl;
    surfaceScalarField alphaPhi
    (
        IOobject
        (
            "alphaPhi",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar(phi.dimensions(), Zero)
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

// Relative Velocity
autoPtr<relativeVelocityModel> UdmModelPtr
(
    relativeVelocityModel::New
    (
        transportProperties,
        mixture,
        "pore:"
    )
);
autoPtr<relativeVelocityModel> fUdmModelPtr 
(
    relativeVelocityModel::New
    (
        transportProperties,
        mixture,
        "free:"
    )
);
Info<<"Done makeing realtive velocity models"<<endl;






#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"

label pRefCell = 0;
scalar pRefValue = 0.0;

if (p.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
}

mesh.setFluxRequired(p.name());
mesh.setFluxRequired(alpha1.name());


// MULES Correction
IOobject alphaPhi10Header
(
    IOobject::groupName("alphaPhi0", alpha1.group()),
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::NO_WRITE
);

// MULES flux from previous time-step
surfaceScalarField alphaPhi10
(
    alphaPhi10Header,
    phi*fvc::interpolate(alpha1)
);




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


Info << "Reading or preparing overpotential field eta \n" << endl;
volScalarField eta
(
	IOobject
	(
		"eta",
		runTime.timeName(),
		mesh,
		IOobject::READ_IF_PRESENT,
		IOobject::NO_WRITE
	),
	mesh,
	dimensionSet(1,2,-3,0,0,-1,0)
);
Info << "Defining currentDistribution INe \n" << endl;
volVectorField INe
(
	IOobject
	(
		"INe",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	fvc::grad(UNe)*poroM.sigma_s_eff()
);

Info << "Defining currentDistribution IPe \n" << endl;
volVectorField IPe
(
	IOobject
	(
		"IPe",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	fvc::grad(UPe)*poroM.sigma_s_eff()
);
Info << "Defining currentDistribution Ie \n" << endl;
volVectorField Ie
(
	IOobject
	(
		"Ie",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	fvc::grad(UPe)*poroM.sigma_s_eff()
);



#include "createMRF.H"
#include "createFvOptions.H"

// massAndSpeciesTransferModel
autoPtr<massAndSpeciesTransferModel> mSTaPtr
(
    massAndSpeciesTransferModel::New
    (
        transportProperties,
        mesh,
        mixture
    )
);

Info<< "Reading speciesTransport\n" << endl;
// speciesTransport
speciesTransport sTp(mesh,transportProperties,mSTaPtr);

// surface coverage
autoPtr<coverageModel> thetaModelPtr
(
    coverageModel::New
    (
        transportProperties,
        mesh,
        mixture
    )
);
// concentrationloss
autoPtr<concentrationLossModel> CLModelPtr
(
    concentrationLossModel::New
    (
        transportProperties,
        mesh,
        mSTaPtr
    )
);
//mixture.correct(mSTaPtr->p_water());
Info<<"Creating Darcy drag field"<<endl;
//Info<<min(poroM.K())<<endl;
//Info<<max(poroM.K())<<endl;
volScalarField Drd
(
        IOobject
        (
            "Drd",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimMass/dimTime/dimLength/dimLength/dimLength,1)
);

    //Labeling the places that have porous media
    forAll(poroM.eps(),celli)
    {
        if(poroM.eps()[celli]>0.99)
        {
            Solid[celli]=0;
        }
        else
        {
            Solid[celli]=1;
        }
    }
    Solid=max(Solid-poroM.PeC()-poroM.NeC(),0.0);
    forAll(mesh.boundary(), patchi)
    {
        forAll(Solid.boundaryField()[patchi],facei)
        {
            if(poroM.eps().boundaryField()[patchi][facei]>0.99 or (poroM.PeC()+poroM.NeC())().boundaryField()[patchi][facei]>=1.0)
            {
                Solid.boundaryFieldRef()[patchi][facei]=0;
            }
            else
            {
                Solid.boundaryFieldRef()[patchi][facei]=1;
            }
        }
    }
surfaceScalarField phir(fvc::flux((1-Solid)*fUdmModelPtr->Udm()+Solid*UdmModelPtr->Udm()));

//create constant fields for rdrPiso and correction of phir
volScalarField chi1
(
        IOobject
        (
            "chi1",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless,1)
);
surfaceScalarField theta1
(
        IOobject
        (
            "theta1",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless,1)
);
surfaceScalarField w1
(
        IOobject
        (
            "w1",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        1.0/Foam::mag(mesh.delta())
);
//starts with rdrPiso
{
const surfaceScalarField absDEpsf = Foam::mag(fvc::snGrad(poroM.eps())*mesh.magSf());  // φ_nb - φ_P  with correct sign
// owner–neighbour addressing


// final normalisation χₚ = Σ_f |Δφ| / max(Σ_f |Δφ|, ξ)
chi1=fvc::surfaceSum(absDEpsf)*dimensionedScalar(dimless/dimLength,1);
forAll(chi1, cellI)
{
    const scalar sum = chi1[cellI];
    chi1[cellI] = sum/(max(sum,VSMALL));          // VSMALL ≈ 1e-15
}

surfaceScalarField absOmegaf(Foam::mag(fvc::snGrad(chi1*poroM.eps())*mesh.magSf()));  // Δ(Δφ)_f  := dPhi_nb − dPhi_P

forAll(absOmegaf, fI)
{
    theta1[fI] = absOmegaf[fI]/(max(absOmegaf[fI], VSMALL));
}

}
surfaceScalarField wgtf(w1*theta1);
volScalarField denP(max(fvc::surfaceSum(wgtf),dimensionedScalar(dimless/dimLength,VSMALL)));
/*
// prepare for phir or Udm
volScalarField chi2
(
        IOobject
        (
            "chi2",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless,1)
);

surfaceScalarField theta2
(
        IOobject
        (
            "theta2",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless,1)
);


{
const surfaceScalarField absDSolidf = Foam::mag(fvc::snGrad(Solid)*mesh.magSf());  // φ_nb - φ_P  with correct sign
// owner–neighbour addressing

// final normalisation χₚ = Σ_f |Δφ| / max(Σ_f |Δφ|, ξ)
chi2=fvc::surfaceSum(absDSolidf)*dimensionedScalar(dimless/dimLength,1);
forAll(chi2, cellI)
{
    const scalar sum = chi2[cellI];
    chi2[cellI] = sum/(max(sum,VSMALL));          // VSMALL ≈ 1e-15
}

surfaceScalarField absOmegaf(Foam::mag(fvc::snGrad(chi2*Solid)*mesh.magSf()));  // Δ(Δφ)_f  := dPhi_nb − dPhi_P

forAll(absOmegaf, fI)
{
    theta2[fI] = absOmegaf[fI]/(max(absOmegaf[fI], VSMALL));
}
}
//Info<<"before List"<<endl;
static labelList activeCells2;
if (activeCells2.empty())
{
    DynamicList<label> tmp;
    forAll(chi2, cellI)
    {
        if (chi2[cellI] >= 0.5) tmp.append(cellI);
    }
    activeCells2.transfer(tmp);   // shrink-to-fit
    Info<< "Active-cell list size : " << activeCells2.size() << nl << endl;
}*/
static labelList activeCells1;
if (activeCells1.empty())
{
    DynamicList<label> tmp;
    forAll(chi1, cellI)
    {
        if (chi1[cellI] >= 0.5) tmp.append(cellI);
    }
    activeCells1.transfer(tmp);   // shrink-to-fit
    Info<< "Active-cell list size : " << activeCells1.size() << nl << endl;
}
//Info<<"after List"<<endl;


