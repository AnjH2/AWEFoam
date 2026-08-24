{
    phir=fvc::flux(UdmModel.Udm());
    
    

  for (int aCorr=0; aCorr<nAlphaCorr; aCorr++)
  {
        // Gas-phase balance: mixture convection, interphase mass transfer,
        // bubble dispersion and relative (drift) transport.
        fvScalarMatrix alpha1Eqn
        (
            fvm::ddt(epsilon1*rho1,alpha1)
            +fvm::div(phi*fvc::interpolate(rho1),alpha1,"div(U,alphad)")
            ==

            -mSTa.mDot()[1]
            +fvm::laplacian((1-Mem)*epsilon1*rho1*UdmModel.Ddm(), alpha1)
            -fvm::div(phir*fvc::interpolate((1-Mem+VSMALL)*rho1),alpha1,"div(Udm,alphad)")
        );
        alpha1Eqn.solve(mesh.solver("alpha1Eqn"));
	
        alphaPhi = alpha1Eqn.flux()/fvc::interpolate(rho1);
        
        alpha2 = 1.0 - alpha1;
        
        
	
	
        Info<< "Phase-1 volume fraction = "
            << alpha1.weightedAverage(mesh.Vsc()).value()
            << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
            << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
            << endl;
    }
    // Reconstruct mixture mass flux after updating the gas phase fraction.
    rhoPhi = fvc::interpolate(mixture.rho())*(alphaPhi*1/fvc::interpolate(alpha1)-phir);
    rho = mixture.rho();
    
    
}
