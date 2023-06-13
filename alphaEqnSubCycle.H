{
    surfaceScalarField alphaPhi
    (
        IOobject
        (
            "alphaPhi",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar(phi.dimensions(), Zero)
    );
    //Psi_alpha=-1*(Psi_si[0] + Psi_si[1])/mixture.rhod();
   Psi_gas=-1*(
        Psi_si[0]/(MW[0]*mixture.p_num()/(R*T)) 
      + Psi_si[1]/(MW[1]*mixture.p_num()/(R*T))
      	);
    Psi_liquid=-1*(
        Psi_si[0]/(rho2) 
      + Psi_si[1]/(rho2)
      	);
    surfaceScalarField phir(fvc::flux(UdmModel.Udm()));

    if (nAlphaSubCycles > 1)
    {
        dimensionedScalar totalDeltaT = runTime.deltaT();
        surfaceScalarField alphaPhiSum
        (
            IOobject
            (
                "alphaPhiSum",
                runTime.timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar(phi.dimensions(), Zero)
        );

        for
        (
            subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
            !(++alphaSubCycle).end();
        )
        {
            #include "alphaEqn.H"
            alphaPhiSum += (runTime.deltaT()/totalDeltaT)*alphaPhi;
        }

        alphaPhi = alphaPhiSum;
    }
    else
    {
        #include "alphaEqn.H"
    }

    // Apply the diffusion term separately to allow implicit solution
    // and boundedness of the explicit advection
    {
        fvScalarMatrix alpha1Eqn
        (
            (fvm::ddt(epsilon1,alpha1) - fvc::ddt(epsilon1,alpha1))
          - fvm::laplacian(turbulence->nut(), alpha1)
          
        );
	alpha1Eqn += Psi_gas;
	
        alpha1Eqn.solve(mesh.solver("alpha1Diffusion"));

        alphaPhi += alpha1Eqn.flux();
        alpha2 = 1.0 - alpha1;

        Info<< "Phase-1 volume fraction = "
            << alpha1.weightedAverage(mesh.Vsc()).value()
            << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
            << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
            << endl;
    }

    rhoPhi = alphaPhi*(fvc::interpolate(rho1) - fvc::interpolate(rho2)) + phi*fvc::interpolate(rho2);
    rho = mixture.rho();
    alpha1Phi=alphaPhi;
    
}
