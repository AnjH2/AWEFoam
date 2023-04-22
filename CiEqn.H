// solve the transport eqaution for species
dimensionedScalar t_OH=u_OH/(u_K+u_OH);
//iterates over each specoes C_[i]
for(label i=0; i<C.size(); i++)
{
	
	volScalarField & Ci=C[i];
	volScalarField & si=s[i];
	volScalarField & Si=S[i];
	
	//dimensionedScalar & 
	//checks if tis a charged species
	if(z[i].value() == 0) //not charged species C_[i]
	{
	Si=-(si/2)*(J)/F;
	//n=2
	}
	else //charged species C_[i]
	{
	Si=-(t_OH/z[i].value()+si/2)*(J)/F;
	//n=2
	}
	//defining species_i transport equation
	fvScalarMatrix CiEqn
	(
        fvm::ddt(epsilon1,Ci)
        +fvm::div(phi,Ci)
        ==
        fvm::laplacian(D_eff[i],Ci)
        +Si // source term
        );
        //relax equation and solves.
        CiEqn.relax();
        solve(CiEqn);
}










        	
