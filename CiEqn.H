// solve the transport eqaution for species in phase 2
//dimensionedScalar t_OH=u_OH/(u_K+u_OH);***********
//iterates over each species C_[i]


for(label i=0; i<C2.size(); i++)
{

	//dimensionedScalar & 
	//checks if it is a charged species, and calculate electrochemical source term.
	if(z[i].value() == 0) //not charged species C_[i]
	{
	react.Psi()[i]=-((react.psi()[i]*Ne+react.psi()[i+4]*Pe)/2)*(J)/F;
	//n=2
	}
	else //charged species C_[i]
	{
	react.Psi()[i]=-(react.t_OH()/z[i].value()+(react.psi()[i]*Ne+react.psi()[i+4]*Pe)/2)*(J)/F;
	//n=2
	}
	//
	fvScalarMatrix C2Eqn
	(
        	fvm::ddt(MW[i]*alpha2*epsilon1,C2[i])
        	+fvm::div(alpha2Phi*MW[i],C2[i])
        	==
        	fvm::laplacian(poroM.eps()*mixture.rho()*alpha2*mixture.D2_eff()[i]*MW[i]/rho2,C2[i])
        	+react.Psi()[i]*MW[i]*poroM.eps()	// electrochemcial source term ***********
        );
	/*if(i<=1)
	{
		fvScalarMatrix Psi_s
		(
		fvm::SuSp(1/react.tau_c()*MW[i],C2[i])-1/react.tau_c()*MW[i]*react.C_sat()[i]
		);
		C2Eqn += Psi_s;
	}*/
	//defining species_i transport equation

        //relax equation and solves.
        C2Eqn.relax();
        solve(C2Eqn);
}

for(label i=0; i<C1.size(); i++)
{

	//defining species_i transport equation
	fvScalarMatrix C1Eqn
	(
        fvm::ddt(MW[i]*alpha1*epsilon1,C1[i])
        +fvm::div(alpha1Phi*MW[i],C1[i])
        //==
        //fvm::laplacian(poroM.eps()*mixture.rho()*alpha1*mixture.D1_eff()*MW[i]/rho1,C1[i])
        //-(fvc::Sp(1/react.tau_c()*MW[i],C2[i])-1/react.tau_c()*MW[i]*react.C_sat()[i]) // source term ***********
        );
        //relax equation and solves.
        C1Eqn.relax();
        solve(C1Eqn);
}








        	
