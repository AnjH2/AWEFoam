    // Solve the Momentum equation

    MRF.correctBoundaryVelocity(U);

    fvVectorMatrix UEqn
    (
        fvm::ddt(mixture.rho(), U) 
      //+ fvm::div(fvc::interpolate(mixture.rho()/pow(epsilon1,2))*phi, U,"div(rhoPhi,U)")
      + fvm::div(fvc::interpolate(1/pow(epsilon1,2))*rhoPhi, U,"div(rhoPhi,U)")
      + fvc::div(
                 1/(pow(epsilon1,2)*mixture.rho())*(alpha1*mixture.rhod())*(alpha2*mixture.rhoc())*Ur*Ur,
                    "div(tauDm)"
                )
      + 1/epsilon1*turbulence->divDevRhoReff(U)
      //- fvm::Sp(fvc::ddt(epsilon1,rho) + (1/pow(epsilon1,2))*fvc::div(rhoPhi), U)
     ==
     fvOptions(mixture.rho(), U)
    );

    UEqn.relax();

    fvOptions.constrain(UEqn);

   /* if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
          ==
            fvc::reconstruct
            (
                (
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                )*mesh.magSf()
            )
        );

        fvOptions.correct(U);
    }*/
