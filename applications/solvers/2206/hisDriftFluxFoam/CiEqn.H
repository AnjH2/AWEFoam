// solve the transport eqaution for species in phase 2
//dimensionedScalar t_OH=u_OH/(u_K+u_OH);***********
//iterates over each species C_[i]
word alphaScheme("div(phi,alpha)");
word alpharScheme("div(phirb,alpha)");

surfaceScalarField alpha2Phi
(
	//phi*fvc::interpolate(alpha2)+fvc::flux(alpha2*((UdmModel.Dcm()&fvc::grad(alpha2))+UdmModel.Ucm()))
	(fvc::flux(phi,rho,alphaScheme) + fvc::flux(-alphaPhi,mixture.rhod(),alpharScheme))/mixture.rhoc()
	   
      //((phi*fvc::interpolate(rho))()-(alphaPhi*fvc::interpolate(mixture.rhod()))())/mixture.rhoc()
        
);

for(label i=0; i<C2.size()-1; i++)
{

	mSTa.correct_Psi_BV(i); //calculates species transfer from exchange current
	//mSTa.Psi_BV()[i].correctBoundaryConditions();
	//

	fvScalarMatrix C2Eqn
	(
        	fvm::ddt(MW[i]*alpha2*epsilon1,C2[i])
        	//+fvm::div(fvc::interpolate(alpha2)*(phi+fvc::flux(UdmModel.Ucm()))*MW[i],C2[i],"div(rhoalpha2*Uc,C2[i])")
        	//+fvm::div(fvc::interpolate(alpha2)*(phi+fvc::flux(UdmModel.Ucm()))*MW[i],C2[i],"div(rhoalpha2*Uc,C2[i])")
        	//+fvm::div(fvc::interpolate(poroM.eps())*alpha2Phi*MW[i],C2[i],"div(rhoalpha2*Uc,C2[i])")
        	+fvm::div(alpha2Phi*MW[i],C2[i],"div(rhoalpha2*Uc,C2[i])")
        	==
        	//-fvm::laplacian(-poroM.eps()*mixture.rho()*alpha2*speciesP.D2_eff()[i]*MW[i]/rho2,C2[i])
        	//+fvm::div(fvc::interpolate(alpha2*UdmModel.Dcm()*MW[i])*fvc::snGrad(alpha2)*mesh.magSf(),C2[i], "div(const*snGrad(alpha2),C2[i])")
        	//+
        	mSTa.Psi_BV()[i]*MW[i]	// electrochemcial source term ***********
        
        );
	
		//calculates concentation of speceies i on the surfaces of the electrodes
	if(i<=1)
	{
		C2Eqn -= -fvm::laplacian(-poroM.eps()*mixture.rho()*alpha2*speciesP.D2_eff()[i]*MW[i]/rho2,C2[i]);
		mSTa.correct_mDot_wall(i,sTp.C2_s(),theta);
		sTp.correct(i,theta);
		Pair<tmp<volScalarField>> mDot = mSTa.mDot(i,1);
    		const volScalarField& mDotc = mDot[0]();
    		const volScalarField& mDotv = mDot[1]();
    		
		//mSTa.correct_Psi_m(i,sTp.C2_s(),theta);
		C2Eqn -= mDotc + mDotv;//mSTa.Psi_m()[i];
		//C2Eqn -= mSTa.ciSource(i);//mSTa.Psi_m()[i];
	} else if (i==2 and mSTa.waterVapour()){
	    C2Eqn -= -fvm::laplacian(-poroM.eps()*mixture.rho()*alpha2*speciesP.D2_eff()[i]*MW[i]/rho2,C2[i]);
		mSTa.correct_mDot_wall(i,sTp.C2_s(),theta);
		sTp.correct(i,theta);
		Pair<tmp<volScalarField>> mDot = mSTa.mDot(i,1);
    		const volScalarField& mDotc = mDot[0]();
    		const volScalarField& mDotv = mDot[1]();
    		
		//mSTa.correct_Psi_m(i,sTp.C2_s(),theta);
		C2Eqn -= mDotc + mDotv;//mSTa.Psi_m()[i];
		
		//C2Eqn -= mSTa.ciSource(i,true);
		//C2Eqn += (mSTa.Psi_m()[0]+mSTa.Psi_m()[1])/(1-mSTa.p_water()/mixture.p_num())-(mSTa.Psi_m()[0]+mSTa.Psi_m()[1]);
	} else if (i>=3){
        Info<<average(speciesP.D2_ambi())<<endl;
		C2Eqn -= -fvm::laplacian(-poroM.eps()*mixture.rho()*alpha2*speciesP.D2_ambi()*MW[i]/rho2,C2[i]);
		C2Eqn -= fvc::div(fvc::flux(speciesP.t()[i]*MW[i]/(z[i]*F)*Ie));
		
	} if (i<=3){
		sTp.correct(i,theta);
		//C2Eqn += -fvm::laplacian(-poroM.eps()*mixture.rho()*alpha2*speciesP.D2_ambi()*MW[i]/rho2,C2[i]);
		//C2Eqn += fvc::div(fvc::flux(speciesP.t_[i]/(speciesP.z_[i]));
		//C2Eqn +=fvc::laplacian(speciesP.D2_eff()[i]*F*z[i]*C2[i]/(R*T),Ue)
		//C2Eqn -=fvm::div((fvc::interpolate(alpha2*speciesP.D2_eff()[i]*F*z[i]*MW[i]/(R*T)))*fvc::snGrad(Ue)*mesh.magSf(),C2[i], "div(const*snGrad(Ue),C2[i])");
		//added poroM.eps()*!!! 06-03-2024
		//romved poroM.eps() its already present from D2_eff() 11-09-2024
	}
    if (i<4){
        C2Eqn.relax();
    
        solve(C2Eqn);
        }
	//Info<<"test3"<<endl;
	
}
//speciesP.correct();
//Info<<average(react.dlnGOH_dlnmKOH(T,react.mKOH(T,C2[3]/1000)))<<endl;
//Info<<average(speciesP.C2_T()/speciesP.C2_0()*(1+react.dlnGOH_dlnmKOH(T,react.mKOH(T,C2[3]/1000))))<<endl;







        	
