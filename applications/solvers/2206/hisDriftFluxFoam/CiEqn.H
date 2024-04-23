// solve the transport eqaution for species in phase 2
//dimensionedScalar t_OH=u_OH/(u_K+u_OH);***********
//iterates over each species C_[i]

surfaceScalarField alpha2Phi
(
	(fvc::interpolate(mixture.rho())*phi-(phi+phir)*fvc::interpolate(alpha1*rho1))/(rho2)
);

for(label i=0; i<C2.size(); i++)
{

	mSTa.correct_Psi_BV(i); //calculates species transfer from exchange current
	//mSTa.Psi_BV()[i].correctBoundaryConditions();
	//

	fvScalarMatrix C2Eqn
	(
        	fvm::ddt(MW[i]*alpha2*epsilon1,C2[i])
        	//+fvm::div(fvc::interpolate(alpha2)*(phi+fvc::flux(UdmModel.Ucm()))*MW[i],C2[i],"div(rhoalpha2*Uc,C2[i])")
        	//+fvm::div(fvc::interpolate(alpha2)*(phi+fvc::flux(UdmModel.Ucm()))*MW[i],C2[i],"div(rhoalpha2*Uc,C2[i])")
        	+fvm::div(alpha2Phi*MW[i],C2[i],"div(rhoalpha2*Uc,C2[i])")
        	==
        	-fvm::laplacian(-poroM.eps()*mixture.rho()*alpha2*speciesP.D2_eff()[i]*MW[i]/rho2,C2[i])
        	//+fvm::div(fvc::interpolate(alpha2*UdmModel.Dcm()*MW[i])*fvc::snGrad(alpha2)*mesh.magSf(),C2[i], "div(const*snGrad(alpha2),C2[i])")
        	+mSTa.Psi_BV()[i]*MW[i]	// electrochemcial source term ***********
        
        );
	
	sTp.correct(i,theta);	//calculates concentation of speceies i on the surfaces of the electrodes
	if(i<=1)
	{
		mSTa.correct_mDot_wall(i,sTp.C2_s(),theta);
		Pair<tmp<volScalarField>> mDot = mSTa.mDot(i,1);
    		const volScalarField& mDotc = mDot[0]();
    		const volScalarField& mDotv = mDot[1]();
    		
		//mSTa.correct_Psi_m(i,sTp.C2_s(),theta);
		C2Eqn -= mDotc + mDotv;//mSTa.Psi_m()[i];
		//C2Eqn -= mSTa.ciSource(i);//mSTa.Psi_m()[i];
	} else if (i==2 and mSTa.waterVapour()){
		mSTa.correct_mDot_wall(i,sTp.C2_s(),theta);
		
		Pair<tmp<volScalarField>> mDot = mSTa.mDot(i,1);
    		const volScalarField& mDotc = mDot[0]();
    		const volScalarField& mDotv = mDot[1]();
    		
		//mSTa.correct_Psi_m(i,sTp.C2_s(),theta);
		C2Eqn -= mDotc + mDotv;//mSTa.Psi_m()[i];
		
		//C2Eqn -= mSTa.ciSource(i,true);
		//C2Eqn += (mSTa.Psi_m()[0]+mSTa.Psi_m()[1])/(1-mSTa.p_water()/mixture.p_num())-(mSTa.Psi_m()[0]+mSTa.Psi_m()[1]);
	} //if (i<=3){
		//C2Eqn +=fvc::laplacian(speciesP.D2_eff()[i]*F*z[i]*C2[i]/(R*T),Ue)
		//C2Eqn -=fvm::div((fvc::interpolate(poroM.eps()*speciesP.D2_eff()[i]*F*z[i]*MW[i]/(R*T)))*fvc::snGrad(Ue)*mesh.magSf(),C2[i], "div(const*snGrad(Ue),C2[i])");//added poroM.eps()*!!! 06-03-2024
	
	//}
        C2Eqn.relax();

        solve(C2Eqn);
	//Info<<"test3"<<endl;
}









        	
