// Dissolved-species transport in the liquid phase (phase 2).
// Couples convection/diffusion to electrochemical production and gas-liquid
// mass transfer for H2, O2 and optional water vapour.

surfaceScalarField alpha2Phi
(
	phi*fvc::interpolate(alpha2)+fvc::flux(alpha2*((UdmModel.Dcm()&fvc::grad(alpha2))+UdmModel.Ucm()))
);

for(label i=0; i<C2.size(); i++)
{
    // Convert local transfer current into the molar source of species i.
	mSTa.correct_Psi_BV(i); 
	
	fvScalarMatrix C2Eqn
	(
        	fvm::ddt(MW[i]*alpha2*epsilon1,C2[i])
        	
        	+fvm::div(alpha2Phi*MW[i],C2[i],"div(rhoalpha2*Uc,C2[i])")
        	==
        	-fvm::laplacian(-poroM.eps()*mixture.rho()*alpha2*speciesP.D2_eff()[i]*MW[i]/rho2,C2[i])
        	
        	+mSTa.Psi_BV()[i]*MW[i]	// electrochemcial source term ***********
        
        );
	
	// H2 and O2: update attached-bubble transfer and subtract the resulting
	// gas-liquid mass-transfer source from the dissolved-species equation.
	if(i<=1)
	{
		
		mSTa.correct_mDot_wall(i,sTp.C2_s(),theta);
		sTp.correct(i,theta);
		Pair<tmp<volScalarField>> mDot = mSTa.mDot(i,1);
    		const volScalarField& mDotc = mDot[0]();
    		const volScalarField& mDotv = mDot[1]();
    		
		
		C2Eqn -= mDotc + mDotv;//mSTa.Psi_m()[i];
		
	} else if (i==2 and mSTa.waterVapour()){
	    // Water transfer is tied to saturation of the generated gas with vapour.
		mSTa.correct_mDot_wall(i,sTp.C2_s(),theta);
		sTp.correct(i,theta);
		Pair<tmp<volScalarField>> mDot = mSTa.mDot(i,1);
    		const volScalarField& mDotc = mDot[0]();
    		const volScalarField& mDotv = mDot[1]();
    		
		
		C2Eqn -= mDotc + mDotv;//mSTa.Psi_m()[i];
		
		
	} if (i<=3){
	    // Update electrode-surface concentrations used by the reaction model.
		sTp.correct(i,theta);
		
	}
        C2Eqn.relax();

        solve(C2Eqn);
}









        	
