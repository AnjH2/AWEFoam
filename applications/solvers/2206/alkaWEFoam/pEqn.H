// Pressure correction for the porous two-phase mixture.
// Enforces continuity and updates the conservative mixture flux and velocity.
{
	
    volScalarField rAU("rAU", 1.0/UEqn.A());

    surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));

    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh));
    // Buoyancy contribution written as a face flux.
    surfaceScalarField phig(-rhorAUf*ghf*fvc::snGrad(rho)*mesh.magSf());

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        (
            fvc::flux(rho*HbyA)
            
        
        )
        + phig
    );

    

    constrainPressure(p_rgh, rho, U, phiHbyA, rhorAUf, MRF);

    
    
    while (pimple.correctNonOrthogonal())
    {
        // Pressure equation includes density changes caused by the evolving
        // phase distribution in the porous volume.
        fvScalarMatrix p_rghEqn
        (
            fvm::laplacian(rhorAUf, p_rgh) 
            == 
            fvc::div(phiHbyA) +fvc::ddt(epsilon1,mixture.rho())
            			
        );
    
        p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));

        p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            // Calculate the conservative fluxes

            phi = (phiHbyA - p_rghEqn.flux())/fvc::interpolate(rho);

            // Explicitly relax pressure for momentum corrector
            p_rgh.relax();

            // Correct the momentum source with the pressure gradient flux
            // calculated from the relaxed pressure

            U = HbyA + rAU*fvc::reconstruct((phig - p_rghEqn.flux())/rhorAUf);
            U.correctBoundaryConditions();
            fvOptions.correct(U);
            
        }
    }
	#include "continuityErrs.H"
    p = p_rgh + rho*gh;

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
}
