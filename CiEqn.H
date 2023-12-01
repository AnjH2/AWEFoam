// solve the transport eqaution for species in phase 2
//dimensionedScalar t_OH=u_OH/(u_K+u_OH);***********
//iterates over each species C_[i]

volScalarField Temp1
(
	IOobject
	(
		"Temp",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
	mesh,
	dimensionSet(0,-3,0,0,1,0,0)
);
/*volScalarField Zero_
(
	IOobject
	(
		"Zero",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
	mesh,
	dimensionSet(0,-3,0,0,1,0,0)
);*/

for(label i=0; i<C2.size(); i++)
{

	mSTa.correct_Psi_BV(i); //calculates species transfer from exchange current
	//mSTa.Psi_BV()[i].correctBoundaryConditions();
	//

	fvScalarMatrix C2Eqn
	(
        	fvm::ddt(MW[i]*alpha2*poroM.eps(),C2[i])
        	+fvm::div((rhoPhi*fvc::interpolate(alpha2/rho)-fvc::flux(UdmModel.Udm()))*MW[i],C2[i])
        	//- fvc::laplacian(turbulence->nut()*pow(epsilon1,poroM.tau())*MW[i]*C2[i],alpha1)
        	==
        	fvm::laplacian(poroM.eps()*rho*alpha2*speciesP.D2_eff()[i]*MW[i]/rho2,C2[i])
        	+mSTa.Psi_BV()[i]*MW[i]*poroM.eps()	// electrochemcial source term ***********
        
        );

	sTp.correct(i,mSTa.Psi_BV(),theta);	//calculates concentation of speceies i on the surfaces of the electrodes
	if(i<=1)
	{
	Temp1=(C2[i]-k_H[i]*mixture.p_num());
	
	Psi_si[i]=	(Pe+Ne)*(
			//(Temp1>zeroField())* // from liquid to gas
				max(1/react.tau_c()*MW[i]*(C2[i]-k_H[i]*mixture.p_num()),1/react.tau_c()*MW[i]*(C2[i]-k_H[i]*mixture.p_num())*0) //(1-alpha1)*epsilon1*
			);


		C2Eqn += Psi_si[i];
	} 
	//defining species_i transport equation

        //relax equation and solves.

        C2Eqn.relax();

        solve(C2Eqn);
	//Info<<"test3"<<endl;
}

for(label i=0; i<C1.size(); i++)
{
	//Temp1=(C2[i]-k_H[i]*mixture.p_num()[i]);

	//defining species_i transport equation
	fvScalarMatrix C1Eqn
	(
        fvm::ddt(MW[i]*alpha1*epsilon1,C1[i])
        +fvm::div((rhoPhi*fvc::interpolate(alpha1/rho)+fvc::flux(UdmModel.Udm()))*MW[i],C1[i])  
        //- fvc::laplacian(turbulence->nut()*pow(epsilon1,poroM.tau())*MW[i]*C1[i],alpha1)  
        );
        //relax equation and solves.
        C1Eqn -= mSTa.Psi_m()[i];
        C1Eqn.relax();
        solve(C1Eqn);
        }
}








        	
