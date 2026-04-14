{
	
    volScalarField rAUd("rAUd", 1.0/(UEqn.A() + Drd));   // effective inverse diagonal
    volScalarField B("B", UEqn.A()/(UEqn.A() + Drd));                    // = A/(A+Drd)

    surfaceScalarField Bf(fvc::interpolate(B));            // <-- consistent face B
    surfaceScalarField rhorAUdf = fvc::interpolate(mixture.rho()*rAUd, "interpolate(rhorAUd)");

    volVectorField HbyA(constrainHbyA(rAUd*UEqn.H(), U, p_rgh));
    surfaceScalarField  HbyAf(fvc::flux(HbyA));


    // (Option 1) simplest: B already “inside” HbyA via rAUd, so don’t multiply by Bf again:
    surfaceScalarField phiHbyA("phiHbyA", fvc::flux(mixture.rho()*HbyA));
    tmp<volScalarField> rAtUd(rAUd);

    if (pimple.consistent())
    {
        rAtUd = 1.0/max(1.0/rAUd - UEqn.H1(), 0.1/rAUd);

        // IMPORTANT: use p_rgh, not p
        // and keep your buoyancy treatment (phig) unchanged.
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
    //fvc::makeRelative(phiHbyA, rho, U);

    /*while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvc::div(phiHbyA)
          - fvm::laplacian(rAUf, p)
          ==
            (vDotc+vDotv) - fvc::div(fvc::flux(alpha2*UdmModel.Ucm()+alpha1*UdmModel.Udm()))
        );

        pEqn.setReference(pRefCell, pRefValue);

        pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA + pEqn.flux();

            U = HbyA + rAU*fvc::reconstruct((phig + pEqn.flux())/rAUf);
            U.correctBoundaryConditions();
            fvOptions.correct(U);
        }
    }*/
    /*
        while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqn
        (
             
            //-fvm::laplacian(Bf*rhorAUf, p_rgh) +fvc::div(phiHbyA) +fvc::ddt(epsilon1,mixture.rho())
            fvc::div(phiHbyA) - fvm::laplacian(Bf*rhorAUf, p_rgh)
            //+fvc::ddt(epsilon1,mixture.rho())
            +fvc::ddt(rho1-rho2,alpha1)
            			
        );*/
        
    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqn
        (
             
            //-fvm::laplacian(Bf*rhorAUf, p_rgh) +fvc::div(phiHbyA) +fvc::ddt(epsilon1,mixture.rho())
            fvc::div(phiHbyA)
            - fvm::laplacian(rhorAtUdf, p_rgh)== 
            -fvc::ddt(epsilon1,mixture.rho())
            //+fvc::ddt(rho1-rho2,alpha1)
            	
            	
        );
    
        p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));

        p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            // Calculate the conservative fluxes

            rhoPhi = phiHbyA + p_rghEqn.flux();
            //phi = rhoPhi / fvc::interpolate(mixture.rho());
            phi = (phiHbyA + p_rghEqn.flux())/ fvc::interpolate(mixture.rho());
            p_rgh.relax();

            //U = HbyA + rAUd * fvc::reconstruct((phig + p_rghEqn.flux())/rhorAUdf);
            U = fvc::reconstruct(phi)*chi1 +(1-chi1)*(HbyA + rAtUd()*fvc::reconstruct((phig + p_rghEqn.flux())/rhorAtUdf));
            //U=fvc::reconstruct(phi);
            U.correctBoundaryConditions();
            fvOptions.correct(U);

            
        }
    }
    //U=fvc::reconstruct(phi);
	#include "continuityErrs.H"
    p = p_rgh + mixture.rho()*gh;
    //p = -min(p,dimensionedScalar(dimPressure,0))+p;
    //p_rgh=p - rho*gh;
    
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
