    // Solve the Momentum equation

    MRF.correctBoundaryVelocity(U);

    fvVectorMatrix UEqn
    (
        fvm::ddt(mixture.rho(), U) 
      + fvm::div(fvc::interpolate(1/pow(epsilon1,2))*rhoPhi, U,"div(rhoPhi,U)")
      + fvc::div(
                 alpha1*mixture.rhod()*pow(alpha2*mixture.rhoc()/mixture.rho(),2)*Ur*Ur
                +alpha2*mixture.rhoc()*pow(alpha1*mixture.rhod()/mixture.rho(),2)*Ur*Ur, 
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
