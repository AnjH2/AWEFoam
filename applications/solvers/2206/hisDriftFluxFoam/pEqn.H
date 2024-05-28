{
	
    volScalarField rAU("rAU", 1.0/UEqn.A());

    surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));

    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh));

    surfaceScalarField phig(-rhorAUf*ghf*fvc::snGrad(rho)*mesh.magSf());

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        (
            fvc::flux(rho*HbyA)
            //+ MRF.zeroFilter(rhorAUf*fvc::ddtCorr(rho, U, phi))
        
        )
        + phig
    );

    //MRF.makeRelative(fvc::interpolate(rho), phiHbyA);

    constrainPressure(p_rgh, rho, U, phiHbyA, rhorAUf, MRF);

    //fvc::makeRelative(phiHbyA, rho, U);

    /*while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqn
        (
            fvc::div(phiHbyA)
          - fvm::laplacian(rAUf, p_rgh)
          ==
            (vDotc+vDotv) - fvc::div(fvc::flux(alpha2*UdmModel.Ucm()+alpha1*UdmModel.Udm()))
        );

        p_rghEqn.setReference(pRefCell, pRefValue);

        p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA + p_rghEqn.flux();

            U = HbyA + rAU*fvc::reconstruct((phig + p_rghEqn.flux())/rAUf);
            U.correctBoundaryConditions();
            fvOptions.correct(U);
        }
    }*/
    
        while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqn
        (
            fvm::laplacian(rhorAUf, p_rgh) 
            == 
            fvc::div(phiHbyA) +fvc::ddt(mixture.rho())
            			
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
