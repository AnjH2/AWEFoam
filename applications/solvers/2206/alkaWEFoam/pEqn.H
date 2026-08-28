// Pressure correction enforcing mixture continuity (paper Eq. 3). The Darcy
// drag modifies the inverse momentum diagonal, after which conservative mass
// flux and mixture velocity are reconstructed from the pressure correction.
{
	
	// Incorporate Darcy drag in the effective inverse momentum diagonal.
    volScalarField rAUd("rAUd", 1.0/(UEqn.A() + Drd));   

    surfaceScalarField rhorAUdf = fvc::interpolate(mixture.rho()*rAUd, "interpolate(rhorAUd)");

    volVectorField HbyA(constrainHbyA(rAUd*UEqn.H(), U, p_rgh));
    
    surfaceScalarField phiHbyA("phiHbyA", fvc::flux(mixture.rho()*HbyA));
    tmp<volScalarField> rAtUd(rAUd);

    if (pimple.consistent())
    {
        rAtUd = 1.0/max(1.0/rAUd - UEqn.H1(), 0.1/rAUd);

     
        phiHbyA += fvc::interpolate(mixture.rho()*(rAtUd() - rAUd))
                * fvc::snGrad(p_rgh)*mesh.magSf();

        HbyA -= (rAUd - rAtUd())*fvc::grad(p_rgh);
    }
    
    constrainPressure(p_rgh, mixture.rho(), U, phiHbyA, rhorAUdf, MRF);
    surfaceScalarField rhorAtUdf = fvc::interpolate(mixture.rho()*rAtUd(), "interpolate(rhorAtUd)");
    surfaceScalarField phig
    (
        - rhorAtUdf * ghf * fvc::snGrad(mixture.rho()) * mesh.magSf()
    );
    
    phiHbyA += phig;
            
    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqn
        (
             
            
            fvc::div(phiHbyA)
            - fvm::laplacian(rhorAtUdf, p_rgh)== 
            -fvc::ddt(epsilon1,mixture.rho())
            	
            	
        );
    
        p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));

        p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            // Calculate the conservative fluxes

            rhoPhi = phiHbyA + p_rghEqn.flux();
            phi = (phiHbyA + p_rghEqn.flux())/ fvc::interpolate(mixture.rho());
            p_rgh.relax();

            //U = HbyA + rAUd * fvc::reconstruct((phig + p_rghEqn.flux())/rhorAUdf);
            // At porosity transitions chi1 selects direct flux reconstruction;
            // elsewhere use the pressure-corrected momentum predictor.
            U = fvc::reconstruct(phi)*chi1 +(1-chi1)*(HbyA + rAtUd()*fvc::reconstruct((phig  + p_rghEqn.flux())/rhorAtUdf));
            U.correctBoundaryConditions();
            fvOptions.correct(U);

            
        }
    }
	#include "continuityErrs.H"
    p = p_rgh + mixture.rho()*gh;
    
    if (p_rgh.needReference())
    {
        p_rgh += dimensionedScalar
        (
            "p_rgh",
            p_rgh.dimensions(),
            pRefValue - getRefCellValue(p_rgh, pRefCell)
        );
    }
}
