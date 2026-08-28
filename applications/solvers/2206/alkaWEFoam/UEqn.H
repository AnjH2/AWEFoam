// Porous-medium mixture momentum equation (paper Eq. 5). The matrix contains
// transient and convective momentum, drift stress and viscous/turbulent
// stresses; Darcy resistance is coupled through the pressure correction.
    // Solve the Momentum equation

    MRF.correctBoundaryVelocity(U);

    // Assemble the mixture momentum operator. The explicit Ur*Ur term is
    // the drift-induced stress from unresolved phase-relative motion.
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
      
     ==
     fvOptions(mixture.rho(), U)
    );

    UEqn.relax();

    fvOptions.constrain(UEqn);

