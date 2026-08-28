// Liquid-phase species transport (paper Eq. 8). The solved species are
// H2, O2, H2O and OH-; K+ is recovered from electroneutrality. The
// equations couple advection, effective diffusion, Butler-Volmer production
// and gas-liquid mass transfer.


// C2.size()-1 deliberately skips K+: the binary electrolyte is closed by
// electroneutrality, so only OH- is transported explicitly.
for(label i=0; i<C2.size()-1; i++)
{
	// Convert the local transfer current into the molar source for species i
	// using reaction stoichiometry and Faraday's law (paper Eq. 13).
	mSTa.correct_Psi_BV(i); //calculates species transfer from exchange current
	
	fvScalarMatrix C2Eqn
	(
        	fvm::ddt(MW[i]*alpha2*epsilon1,C2[i])
        	+fvm::div((rhoPhi-fvc::interpolate(rho1)*alphaPhi)/rho2*MW[i],C2[i],"div(rhoalpha2*Uc,C2[i])")
        	==
        	mSTa.Psi_BV()[i]*MW[i]	// electrochemcial source term ***********
        
        );
	
	// Neutral H2/O2 use molecular diffusion and exchange mass with attached
	// and detached bubbles through the selected mass-transfer model.
	if(i<=1)
	{
		C2Eqn -= -fvm::laplacian(-poroM.eps()*mixture.rho()*alpha2*speciesP.D2_eff()[i]*MW[i]/rho2,C2[i],"laplacian(D,C[i])");
		mSTa.correct_mDot_wall(i,sTp.C2_s(),theta);
		sTp.correct(i,theta);
		Pair<tmp<volScalarField>> mDot = mSTa.mDot(i,1);
    		const volScalarField& mDotc = mDot[0]();
    		const volScalarField& mDotv = mDot[1]();
    		
		
		C2Eqn -= mDotc + mDotv;
	} else if (i==2 and mSTa.waterVapour()){
	    C2Eqn -= -fvm::laplacian(-poroM.eps()*mixture.rho()*alpha2*speciesP.D2_eff()[i]*MW[i]/rho2,C2[i],"laplacian(D,C[i])");
		mSTa.correct_mDot_wall(i,sTp.C2_s(),theta);
		sTp.correct(i,theta);
		Pair<tmp<volScalarField>> mDot = mSTa.mDot(i,1);
    		const volScalarField& mDotc = mDot[0]();
    		const volScalarField& mDotv = mDot[1]();
    		
		
		C2Eqn -= mDotc + mDotv;
		
		
	} else if (i>=3){
        // OH- uses the effective ambipolar diffusivity. The current-divergence
        // term represents ionic migration consistently with the current field.
		C2Eqn -= -fvm::laplacian(-poroM.eps()*mixture.rho()*alpha2*speciesP.D2_ambi()*MW[i]/rho2,C2[i],"laplacian(D,C[i])");
		C2Eqn -= fvc::div(fvc::flux(MW[i]/(z[i]*F)*Ie));
		
		
	} if (i<=3){
		sTp.correct(i,theta);
		
	}
    if (i<4){
        C2Eqn.relax();
    
        solve(C2Eqn);
        }
	
	
}







        	
