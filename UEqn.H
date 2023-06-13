    // Solve the Momentum equation

    MRF.correctBoundaryVelocity(U);

    fvVectorMatrix UEqn
    (
        fvm::ddt(epsilon1*rho, U) + fvm::div(rhoPhi, U)
      + MRF.DDt(epsilon1*rho, U)
      + fvc::div(UdmModel.tauDm())
      + turbulence->divDevRhoReff(U)
     ==
        fvOptions(epsilon1*rho, U)
    );

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (pimple.momentumPredictor())
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
    }
