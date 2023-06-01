// solve the current distribution in both electrolyte and electrode.
//rewritting butler volmer eqaution so Ue, UNe or UPe can be solved fully implicit.
		//index refers to the given species.
		//[0]->H2; [2]->H2O; [1]->O2; [3]->OH
//******************************************* constants used in eqation
		volScalarField C1Ne
		(
		Ne*as[0]*react.i0()[0]*pow(1-alpha1,n)*Foam::pow(C2[3]/C2_ref[3],2)*(C2[0]/C2_ref[0])
		);
		volScalarField C1Pe
		(
		Pe*as[1]*react.i0()[1]*pow(1-alpha1,n)*Foam::pow(C2[3]/C2_ref[3],2)
		);
		volScalarField C2Ne
		(
		Ne*as[0]*react.i0()[0]*pow(1-alpha1,n)*Foam::pow(C2[2]/C2_ref[2],2)
		);
		volScalarField C2Pe
		(
		Pe*as[1]*react.i0()[1]*pow(1-alpha1,n)*Foam::pow(C2[1]/C2_ref[1],0.5)*(C2[2]/C2_ref[2])
		);
		//opdated 1511-2022
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

        	fvScalarMatrix UeEqn
        	(fvm::laplacian((Foam::pow(z[3]*F,2)*mixture.D2_eff()[3]*C2[3])/(R*T), Ue)+fvc::laplacian(z[3]*F*mixture.D2_eff()[3],C2[3])	
        	==-1*(Je_Pe+Je_Ne));//fvm::Sp(-1,JPe)
        	
        	//define potential equation for negative electrode.
        	 //region 2
        	fvScalarMatrix UNeEqn
        	(fvm::laplacian(max(poroM.sigma_s_eff()*Ne,poroM.sigma_s_eff()*1e-10), UNe)	==J_Ne	);
        	//define potential equation for negative electrode.
        	//region 4
        	fvScalarMatrix UPeEqn
        	(fvm::laplacian(max(poroM.sigma_s_eff()*Pe,poroM.sigma_s_eff()*1e-10), UPe)	==J_Pe	);
        	/*//relax equations, might be unnecessary as they are not time dependt.
        	UNeEqn.relax();
        	UeEqn.relax();
        	UPeEqn.relax();
		*/
        	//solve eqations.
        	{
            	UNeEqn.solve();
            	UeEqn.solve();
            	UPeEqn.solve();
            	}
            	//evaluates J for the solving of C_[i] and plotting.

            	eta=(UNe-Ue-react.E0()[0])*Ne+(UPe-Ue-react.E0()[1])*Pe;
        	
        	J=
        	Pe*as[1]*react.i0()[1]*pow(1-alpha1,n)*(Foam::pow(C2[3]/C2_ref[3],2)*exp((react.alphaT()[2]*F*eta)/(R*T))-Foam::pow(C2[1]/C2_ref[1],0.5)*(C2[2]/C2_ref[2])*exp((-react.alphaT()[3]*F*eta)/(R*T)))+
        	Ne*as[0]*react.i0()[0]*pow(1-alpha1,n)*(Foam::pow(C2[3]/C2_ref[3],2)*(C2[0]/C2_ref[0])*exp((react.alphaT()[0]*F*eta)/(R*T))-Foam::pow(C2[2]/C2_ref[2],2)*exp((-react.alphaT()[1]*F*eta)/(R*T)))
    		;
            	
