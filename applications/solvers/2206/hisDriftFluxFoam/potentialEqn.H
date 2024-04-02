// solve the current distribution in both electrolyte and electrode.
//rewritting butler volmer eqaution so Ue, UNe or UPe can be solved fully implicit.
		//index refers to the given species.
		//[0]->H2; [2]->H2O; [1]->O2; [3]->OH
//******************************************* constants used in eqation
		volScalarField C1Ne
		(
		Ne*as[0]*react.i0()[0]*(1-theta)*(CRa_Ne) //Foam::pow(((react.Psi()[3]/as[0])/(mixture.D2_eff()[3])+C2[3])/C2[3],2)
		);

		volScalarField C1Pe
		(
		Pe*as[1]*react.i0()[1]*(1-theta)*(CRa_Pe) //Foam::pow(((react.Psi()[3]/as[1])/(mixture.D2_eff()[3])+C2[3])/C2[3],2)
		);

		volScalarField C2Ne
		(
		Ne*as[0]*react.i0()[0]*(1-theta)*(CRc_Ne) //Foam::pow(((react.Psi()[2]/as[0])/(mixture.D2_eff()[2])+C2[2])/C2[2],2)
		);

		volScalarField C2Pe
		(
		Pe*as[1]*react.i0()[1]*(1-theta)*(CRc_Pe) //(((react.Psi()[2]/as[1])/(mixture.D2_eff()[2])+C2[2])/C2[2])
		);
		//opdated 1511-2022 fixed error
		//opdated 0307-2023 removed H2 and O2 fraction, added concentration gradient to the surface of the electrodes.
		//opdated 2609-2023 added class to model concentration losses and inserted H2 and O2 again.

//******************************************* constants used in eqation
		
		volScalarField C3Ne
		(
		(react.alphaT()[0]*F)/(R*T)
		);
		volScalarField C3Pe
		(
		(react.alphaT()[2]*F)/(R*T)
		);
		volScalarField C4Ne
		(
		(-react.alphaT()[1]*F)/(R*T)
		);
		volScalarField C4Pe
		(
		(-react.alphaT()[3]*F)/(R*T)
		);
//****************************************** for the electrolyte
		volScalarField I1e_Ne
		(
		C1Ne*exp(C3Ne*(UNe-Ue.oldTime()-react.E0()[0]))
		); 
		volScalarField I1e_Pe
		(
		C1Pe*exp(C3Pe*(UPe-Ue.oldTime()-react.E0()[1]))
		);
		volScalarField I2e_Ne
		(
		C2Ne*exp(C4Ne*(UNe-Ue.oldTime()-react.E0()[0]))
		);
		volScalarField I2e_Pe
		(
		C2Pe*exp(C4Pe*(UPe-Ue.oldTime()-react.E0()[1]))
		);

//******************************************** for the electrodes
		volScalarField I1_Ne
		(
		C1Ne*exp(C3Ne*(UNe.oldTime()-Ue-react.E0()[0]))
		); 
		volScalarField I1_Pe
		(
		C1Pe*exp(C3Pe*(UPe.oldTime()-Ue-react.E0()[1]))
		);
		volScalarField I2_Ne
		(
		C2Ne*exp(C4Ne*(UNe.oldTime()-Ue-react.E0()[0]))
		);
		volScalarField I2_Pe
		(
		C2Pe*exp(C4Pe*(UPe.oldTime()-Ue-react.E0()[1]))
		);

	
//*************************************** for the electrolyte
		fvScalarMatrix	Je_Ne
		(
		I1e_Ne*(1+C3Ne*Ue.oldTime())-fvm::Sp(I1e_Ne*C3Ne,Ue)
		+I2e_Ne*(-1-C4Ne*Ue.oldTime())+fvm::Sp(I2e_Ne*C4Ne,Ue)
		);
		fvScalarMatrix Je_Pe
		(
		I1e_Pe*(1+C3Pe*Ue.oldTime())-fvm::Sp(I1e_Pe*C3Pe,Ue)
		+I2e_Pe*(-1+C4Pe*Ue.oldTime())+fvm::Sp(I2e_Pe*C4Pe,Ue)
		);

//****************************************** for the electrodes
		fvScalarMatrix	J_Ne
		(
		I1_Ne*(1+C3Ne*UNe.oldTime())-fvm::Sp(I1_Ne*C3Ne,UNe)
		+I2_Ne*(-1-C4Ne*UNe.oldTime())+fvm::Sp(I2_Ne*C4Ne,UNe)
		);
		fvScalarMatrix J_Pe
		(
		I1_Pe*(1+C3Pe*UPe.oldTime())-fvm::Sp(I1_Pe*C3Pe,UPe)
		+I2_Pe*(-1-C4Pe*UPe.oldTime())+fvm::Sp(I2_Pe*C4Pe,UPe)
		);

        	//define potential equation for electrolyte.
        	//region 2 3 4
        	/*fvScalarMatrix UeEqn
        	(
        	fvm::laplacian(calKappa(T,C2[3]/1000)*(pow(epsilon1,poroM.tau())*pow(alpha2,1.5)), Ue)//-z[3]*F*fvc::div(rhoPhi/fvc::interpolate(rho),C2[3])	
        	==-1*(Je_Pe+Je_Ne)
        	);*/
		fvScalarMatrix UeEqn
        	(
        	fvm::laplacian((Foam::pow(z[3]*F,2)*speciesP.D2_eff()[3]*C2[3])/(R*T), Ue)//+fvc::laplacian(z[3]*F*speciesP.D2_eff()[3],C2[3])	
        	==-1*(Je_Pe+Je_Ne)
        	);
		
		
        	//define potential equation for negative electrode.
        	 //region 2
        	fvScalarMatrix UNeEqn
        	(fvm::laplacian(poroM.sigma_s_eff(), UNe)	==J_Ne	);
        	//(fvm::laplacian(max(poroM.sigma_s_eff(),poroM.sigma_s_eff()*1e-10), UNe)	==J_Ne	);
        	//define potential equation for negative electrode.
        	//region 4
        	fvScalarMatrix UPeEqn
        	(fvm::laplacian(poroM.sigma_s_eff(), UPe)	==J_Pe	);
        	//relax equations, might be unnecessary as they are not time dependt.
        	//UNeEqn.relax();
        	//UeEqn.relax();
        	//UPeEqn.relax();
		
        	//solve eqations.
        	{
            	UNeEqn.solve();

            	UeEqn.solve();

            	UPeEqn.solve();

            	}
            	//evaluates J for the solving of C_[i] and plotting.
            	eta=(UNe-Ue-react.E0()[0])*Ne+(UPe-Ue-react.E0()[1])*Pe;
            	eta.correctBoundaryConditions();
        	J=
        	Pe*as[1]*react.i0()[1]*(1-theta)*(CRa_Pe*exp((react.alphaT()[2]*F*eta)/(R*T))-(CRc_Pe)*exp((-react.alphaT()[3]*F*eta)/(R*T)))+
        	Ne*as[0]*react.i0()[0]*(1-theta)*(CRa_Ne*exp((react.alphaT()[0]*F*eta)/(R*T))-CRc_Ne*exp((-react.alphaT()[1]*F*eta)/(R*T)))
    		;
		J.correctBoundaryConditions();
