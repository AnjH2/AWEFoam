{
    phir=fvc::flux(UdmModel.Udm());
    
      Pair<tmp<volScalarField>> mDot =
        mSTa.mDot();
    const volScalarField& mDotc = mDot[0]();
    const volScalarField& mDote = mDot[1]();
    
	
  for (int aCorr=0; aCorr<nAlphaCorr; aCorr++)
  {
  	
  	
        fvScalarMatrix alpha1Eqn
        (
            fvm::ddt(epsilon1*rho1,alpha1)
            +fvm::div(phi*fvc::interpolate((1-Mem+VSMALL)*rho1),alpha1,"div(rhod*U,alphad)")
            ==
            -mDote//+mDotc
            +fvm::laplacian((1-Mem+VSMALL)*rho1*UdmModel.Ddm(), alpha1)
            -fvm::div(phir*fvc::interpolate((1-Mem+VSMALL)*rho1),alpha1,"div(rhod*Udm,alphad)")
        );

        alpha1Eqn.solve(mesh.solver("alpha1Eqn"));
	
        alphaPhi = alpha1Eqn.flux()/fvc::interpolate(epsilon1*rho1);
        
        alpha2 = 1.0 - alpha1;
        
        
	
	
        Info<< "Phase-1 volume fraction = "
            << alpha1.weightedAverage(mesh.Vsc()).value()
            << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
            << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
            << endl;
    }
	//alphaPhi=alphaPhi*fvc::interpolate(epsilon1); //to make it superficial
    //rhoPhi = alphaPhi*(fvc::interpolate(mixture.rhod()) - rho2) + phi*rho2; //original
    //rhoPhi = fvc::interpolate(alpha1*mixture.rhod())*(phi+phir)+fvc::interpolate(1-alpha1)*(phi+phiL)*rho2; //V665
    //rhoPhi = alphaPhi*(fvc::interpolate(mixture.rhod()) - rho2)+fvc::interpolate(1-alpha1)*(phi+phiL)*rho2; //V722
    //surfaceScalarField alpha2Phi=fvc::interpolate(1-alpha1)*(fvc::interpolate(1/alpha1)*alphaPhi-phir-phiDdm+phiDcm+phiL);
    //rhoPhi = fvc::interpolate(rho1)*alphaPhi+fvc::interpolate(1-alpha1)*(fvc::interpolate(1/alpha1)*alphaPhi-phir-phiDdm+phiDcm+phiL)*rho2; //V754
    rhoPhi = fvc::interpolate(mixture.rho())*(alphaPhi*fvc::interpolate(1/alpha1)-phir*fvc::interpolate(epsilon1));
    rho = mixture.rho();
    
    //rhoPhi = fvc::interpolate(rho)*(fvc::interpolate(1/alpha1)*alphaPhi-phir);; //V755
    
}
