{
	
    volScalarField rAU("rAU", 1.0/UEqn.A());

    surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));
    volScalarField B(1/(1+Drd*rAU));
    surfaceScalarField Bf(fvc::interpolate(B));
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh));
    surfaceVectorField rhoHbyAf( fvc::interpolate(rho*HbyA));
    
    surfaceScalarField phiPc
    (
        (
            //mixture_->surfaceTensionForce()*(1.-Solidf) //activate free fluid interface model
          fvc::interpolate(mixture.pcCoeff()*Solid)*fvc::snGrad(alpha2)//activate porous media model
         - fvc::interpolate(mixture.pcScale()*Solid)*fvc::snGrad(mixture.pc())
        )*rhorAUf*mesh.magSf()
    );
    
    
    surfaceScalarField phig
    (
        -1*fvc::interpolate(Solid)*rhorAUf*ghf*fvc::snGrad(mixture.rhoS())*mesh.magSf()
        -1*fvc::interpolate(1-Solid)*rhorAUf*ghf*fvc::snGrad(mixture.rho())*mesh.magSf()
        //1*fvc::interpolate(Solid)*rhorAUf*fvc::interpolate(mixture.rhoS())*(g & mesh.Sf())
        //+1*fvc::interpolate(1-Solid)*rhorAUf*fvc::interpolate(mixture.rho())*(g & mesh.Sf())
    );
    
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        (
            Bf * (rhoHbyAf & mesh.Sf())
            //+ MRF.zeroFilter(rhorAUf*fvc::ddtCorr(rho, U, phi))
        
        )
        + Bf*(phig+phiPc)
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
             
            fvc::div(phiHbyA) +fvc::ddt(epsilon1,mixture.rho())
            - fvm::laplacian(Bf*rhorAUf, p_rgh) 
            			
        );
    
        p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));

        p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            // Calculate the conservative fluxes

            phi = (phiHbyA + p_rghEqn.flux())/fvc::interpolate(rho);

            // Explicitly relax pressure for momentum corrector
            p_rgh.relax();

            // Correct the momentum source with the pressure gradient flux
            // calculated from the relaxed pressure

            //U = HbyA + rAU*fvc::reconstruct((phig - p_rghEqn.flux())/rhorAUf);
            U=fvc::reconstruct(phi);
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
