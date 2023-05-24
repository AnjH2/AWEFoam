// solve the transport eqaution for species in phase 2
//dimensionedScalar t_OH=u_OH/(u_K+u_OH);***********
//iterates over each species C_[i]

//**** wanted it to be a function in relectiveVelocityModel, but could not make it work ****
tmp<volScalarField> 	betac(alpha2*rho2);
tmp<volScalarField> 	betad(alpha1*rho1);
volVectorField		Ucm(betad*UdmModel.Udm()/betac);

surfaceScalarField phik2(phi-fvc::flux(Ucm));
surfaceScalarField phik1(phi+fvc::flux(UdmModel.Udm()));
for(label i=0; i<C2.size(); i++)
{
	
	//volScalarField& Ci=C2[i];
	//dimensionedScalar& MWi=MW[i];
	//volScalarField & si=mixture.s()[i];
	//volScalarField & Si=S[i];
	
	//dimensionedScalar & 
	//checks if tis a charged species
	/*if(z[i].value() == 0) //not charged species C_[i]
	{
	Si=-(si/2)*(J)/F;
	//n=2
	}
	else //charged species C_[i]
	{
	Si=-(t_OH/z[i].value()+si/2)*(J)/F;
	//n=2
	}*/
	//defining species_i transport equation
	fvScalarMatrix C2Eqn
	(
        fvm::ddt(MW[i]*alpha2*epsilon1,C2[i])
        +fvm::div(phik2*(fvc::interpolate(alpha2))*MW[i],C2[i])
        ==
        fvm::laplacian(poroM.eps()*mixture.rho()*alpha2*mixture.D2_eff()[i]*MW[i]/rho2,C2[i])
        //+Si // source term ***********
        );
        //relax equation and solves.
        C2Eqn.relax();
        solve(C2Eqn);
}










        	
