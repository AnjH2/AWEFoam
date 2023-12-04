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
        	+fvm::div((rhoPhi*fvc::interpolate(alpha2/rho))*MW[i],C2[i])
        	//+fvm::div((rhoPhi*fvc::interpolate(alpha2/rho)-fvc::flux(UdmModel.Udm()))*MW[i],C2[i])
        	//- fvc::laplacian(turbulence->nut()*pow(epsilon1,poroM.tau())*MW[i]*C2[i],alpha1)
        	==
        	fvm::laplacian(poroM.eps()*rho*alpha2*speciesP.D2_eff()[i]*MW[i]/rho2,C2[i])
        	+mSTa.Psi_BV()[i]*MW[i]*poroM.eps()	// electrochemcial source term ***********
        
        );

	sTp.correct(i,mSTa.Psi_BV(),theta);	//calculates concentation of speceies i on the surfaces of the electrodes
	if(i<=1)
	{

		mSTa.correct_Psi_m(i,sTp.C2_s()); //calculate species transfer based on transfer model
		//mSTa.Psi_m()[i].correctBoundaryConditions();
		C2Eqn += mSTa.Psi_m()[i];
	} else if (i==2 and mSTa.waterVapour()){
		mSTa.correct_Psi_m(i,sTp.C2_s());
		
		C2Eqn += mSTa.Psi_m()[i];
		//C2Eqn += (mSTa.Psi_m()[0]+mSTa.Psi_m()[1])/(1-mSTa.p_water()/mixture.p_num())-(mSTa.Psi_m()[0]+mSTa.Psi_m()[1]);
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
        +fvm::div((rhoPhi*fvc::interpolate(alpha1/rho))*MW[i],C1[i])  
        //+fvm::div((rhoPhi*fvc::interpolate(alpha1/rho)+fvc::flux(UdmModel.Udm()))*MW[i],C1[i])  
        //- fvc::laplacian(turbulence->nut()*pow(epsilon1,poroM.tau())*MW[i]*C1[i],alpha1)  
        );
        //relax equation and solves.
        C1Eqn -= mSTa.Psi_m()[i];
        C1Eqn.relax();
        solve(C1Eqn);
        }
}








        	
