{
    phir=fvc::flux(UdmModel.Udm());
    


    volScalarField divU
    (
	fvc::div(phi)
    );
    
    if (nAlphaSubCycles > 1)
    {
        dimensionedScalar totalDeltaT = runTime.deltaT();
        surfaceScalarField rhoPhiSum
        (
            IOobject
            (
                "rhoPhiSum",
                runTime.timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar(rhoPhi.dimensions(), Zero)
        );

        for
        (
            subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
            !(++alphaSubCycle).end();
        )
        {
            #include "alphaEqn.H"
            rhoPhiSum += (runTime.deltaT()/totalDeltaT)*rhoPhi;
        }

        rhoPhi = rhoPhiSum;
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
            fvm::ddt(epsilon1,alpha1) - fvc::ddt(epsilon1,alpha1)
          - fvm::laplacian(alpha1*UdmModel.Ddm(), alpha1)
        );

        alpha1Eqn.solve(mesh.solver("alpha1Diffusion"));
	
        alphaPhi += alpha1Eqn.flux();
        
        alpha2 = 1.0 - alpha1;
        phir += alpha1Eqn.flux()*fvc::interpolate(1/alpha1);
        
	
	
        Info<< "Phase-1 volume fraction = "
            << alpha1.weightedAverage(mesh.Vsc()).value()
            << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
            << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
            << endl;
    }

    //rhoPhi = alphaPhi*(fvc::interpolate(mixture.rhod()) - rho2) + phi*rho2; //original
    //rhoPhi = fvc::interpolate(alpha1*mixture.rhod())*(phi+phir)+fvc::interpolate(1-alpha1)*(phi+phiL)*rho2; //V665
    //rhoPhi = alphaPhi*(fvc::interpolate(mixture.rhod()) - rho2)+fvc::interpolate(1-alpha1)*(phi+phiL)*rho2; //V722
    //surfaceScalarField alpha2Phi=fvc::interpolate(1-alpha1)*(fvc::interpolate(1/alpha1)*alphaPhi-phir-phiDdm+phiDcm+phiL);
    //rhoPhi = fvc::interpolate(rho1)*alphaPhi+fvc::interpolate(1-alpha1)*(fvc::interpolate(1/alpha1)*alphaPhi-phir-phiDdm+phiDcm+phiL)*rho2; //V754
    rhoPhi = fvc::interpolate(mixture.rho())*(alphaPhi*fvc::interpolate(1/alpha1)-phir);
    rho = mixture.rho();
    
    //rhoPhi = fvc::interpolate(rho)*(fvc::interpolate(1/alpha1)*alphaPhi-phir);; //V755
    
}
