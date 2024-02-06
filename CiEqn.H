// solve the transport eqaution for species in phase 2
//dimensionedScalar t_OH=u_OH/(u_K+u_OH);***********
//iterates over each species C_[i]
for(label i=0; i<C2.size(); i++)
{

	mSTa.correct_Psi_BV(i); //calculates species transfer from exchange current
	//mSTa.Psi_BV()[i].correctBoundaryConditions();
	//

	fvScalarMatrix C2Eqn
	(
        	fvm::ddt(MW[i]*alpha2*poroM.eps(),C2[i])
        	+fvm::div(fvc::interpolate(alpha2)*(phi-fvc::interpolate(alpha1*rho1/rho)*fvc::flux(UdmModel.Udm()))*MW[i],C2[i])
        	//+fvm::div(fvc::interpolate(alpha2)*(phi+fvc::interpolate(alpha2*rho2/rho)*fvc::flux(UdmModel.Udm()))*MW[i],C2[i])
        	//+fvm::div(fvc::interpolate(alpha2)*(phi+fvc::interpolate(alpha2)*fvc::flux(UdmModel.Udm()))*MW[i],C2[i])
        	//- fvc::laplacian(turbulence->nut()*pow(epsilon1,poroM.tau())*MW[i]*C2[i],alpha1)
        	==
        	-fvm::laplacian(-poroM.eps()*rho*alpha2*speciesP.D2_eff()[i]*MW[i]/rho2,C2[i])
        	+mSTa.Psi_BV()[i]*MW[i]	// electrochemcial source term ***********
        
        );
	
	sTp.correct(i,theta);	//calculates concentation of speceies i on the surfaces of the electrodes
	if(i<=1)
	{
		mSTa.correct_mDot_wall(i,sTp.C2_s(),theta);
		
		Pair<tmp<volScalarField>> mDot = mSTa.mDot(i);
    		const volScalarField& mDotc = mDot[0]();
    		const volScalarField& mDotv = mDot[1]();

		//mSTa.correct_Psi_m(i,sTp.C2_s(),theta); //calculate species transfer based on transfer model
		
		//mSTa.Psi_m()[i].correctBoundaryConditions();
		C2Eqn -= mDotc + mDotv;//mSTa.Psi_m()[i];
	} else if (i==2 and mSTa.waterVapour()){
		//mSTa.correct_mDot_wall(i,sTp.C2_s(),theta);
	
		Pair<tmp<volScalarField>> mDot = mSTa.mDot(i);
    		const volScalarField& mDotc = mDot[0]();
    		const volScalarField& mDotv = mDot[1]();
    		
		//mSTa.correct_Psi_m(i,sTp.C2_s(),theta);
		C2Eqn -= mDotc + mDotv;//mSTa.Psi_m()[i];
		//C2Eqn -= mSTa.ciSource(i,true);
		//C2Eqn += (mSTa.Psi_m()[0]+mSTa.Psi_m()[1])/(1-mSTa.p_water()/mixture.p_num())-(mSTa.Psi_m()[0]+mSTa.Psi_m()[1]);
	} else if (i==3){
		//C2Eqn +=fvc::laplacian(speciesP.D2_eff()[i]*F*z[i]*C2[i]/(R*T),Ue)
		C2Eqn -=fvm::div((fvc::interpolate(speciesP.D2_eff()[i]*F*z[i]*MW[i]/(R*T)))*fvc::snGrad(Ue)*mesh.magSf(),C2[i], "div(const*snGrad(Ue),C2[i])");
	
	}
        C2Eqn.relax();

        solve(C2Eqn);
	//Info<<"test3"<<endl;
}

for(label i=0; i<C1.size(); i++)
{

	if (i<=1 or mSTa.waterVapour()){

	//defining species_i transport equation
	fvScalarMatrix C1Eqn
	(
        fvm::ddt(MW[i]*alpha1*epsilon1,C1[i])
        +fvm::div(fvc::interpolate(alpha1)*(phi+fvc::interpolate(alpha2*rho2/rho)*fvc::flux(UdmModel.Udm()))*MW[i],C1[i],"div(alphaPhi,"+C1[i].name()+")")//V3
        //V3+fvm::div(fvc::interpolate(alpha1)*(phi-fvc::interpolate(alpha1*rho1/rho)*fvc::flux(UdmModel.Udm()))*MW[i],C1[i],"div(alphaPhi,"+C1[i].name()+")")  
        //+fvm::div(fvc::interpolate(alpha1)*(phi-fvc::interpolate(alpha1*rho1/rho)*fvc::flux(UdmModel.Udm()))*MW[i],C1[i],"div(alphaPhi,"+C1[i].name()+")")  
        //+fvm::div(fvc::interpolate(alpha1)*(phi-fvc::interpolate(alpha1)*fvc::flux(UdmModel.Udm()))*MW[i],C1[i],"div(alphaPhi,"+C1[i].name()+")")  
        //- fvc::laplacian(turbulence->nut()*pow(epsilon1,poroM.tau())*MW[i]*C1[i],alpha1)  
        );
        //relax equation and solves.
        
        Pair<tmp<volScalarField>> mDot = mSTa.mDot(i);
    	const volScalarField& mDotc = mDot[0]();
    	const volScalarField& mDotv = mDot[1]();
        
        C1Eqn += mDotc + mDotv;//mSTa.Psi_m()[i];
        //C1Eqn += mSTa.ciSource(i, false);
        C1Eqn.relax();
        solve(C1Eqn);
        }
}








        	
