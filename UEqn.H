    // Solve the Momentum equation

    MRF.correctBoundaryVelocity(U);

    fvVectorMatrix UEqn
    (
        fvm::ddt(epsilon1*rho, U) + fvm::div(rhoPhi, U)
      + MRF.DDt(rho, U)
      + fvc::div(UdmModel.tauDm())
      + turbulence->divDevRhoReff(U)
     ==
        fvOptions(rho, U)						//-epsilon1*F*z[3]*C2[3]*fvc::grad(Ue)
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
