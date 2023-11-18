{
    //mixture.correct_rhod();
    //rho=mixture.rho();
    //volScalarField epsRho(epsilon1*rho);
    volScalarField rAU("rAU", 1.0/UEqn.A());
    surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh));

    surfaceScalarField phig(-rhorAUf*ghf*fvc::snGrad(rho)*mesh.magSf());//fvc::interpolate(epsilon1)
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(rho*HbyA)
      + MRF.zeroFilter(rhorAUf*fvc::ddtCorr(rho,U, rhoPhi))
      
    );
    MRF.makeRelative(fvc::interpolate(rho), phiHbyA);
    
    //adjustPhi(phiHbyA, U, p_rgh);
    phiHbyA += phig;
    //adjustPhi(phiHbyA, U, p_rgh);

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p_rgh, rho, U, phiHbyA, rhorAUf, MRF);
    //constrainPressure(p_rgh, U, phiHbyA, rAUf, MRF);//original
    while (pimple.correctNonOrthogonal())
    {

    	fvScalarMatrix p_rghEqn
    	(
		fvm::laplacian(rhorAUf, p_rgh) == fvc::div(phiHbyA) +fvc::ddt(epsilon1,rho) 
    	);
	//p_rghEqn += Psi_si[0]+Psi_si[1];


        p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));
        p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {

            rhoPhi = phiHbyA - p_rghEqn.flux();
            phi=rhoPhi/fvc::interpolate(rho);

            p_rgh.relax();

            U = HbyA + rAU*fvc::reconstruct((phig - p_rghEqn.flux())/rhorAUf);//-
            U.correctBoundaryConditions();
            fvOptions.correct(U);
        }
    }
    #include "continuityErrs.H"
    //#include "compressibleContinuityErrs.H"
    p == p_rgh + rho*gh;

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
