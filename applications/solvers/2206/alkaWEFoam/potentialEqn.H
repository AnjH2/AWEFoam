// Coupled solid/electrolyte current-continuity equations (paper Eqs. 12-25).
// Butler-Volmer kinetics are linearized about the previous potential so UNe,
// UPe and Ue can be solved implicitly before reconstructing eta, J and Ie.

// solve the current distribution in both electrolyte and electrode.
//rewritting butler volmer eqaution so Ue, UNe or UPe can be solved fully implicit.
		//index refers to the given species.
		//[0]->H2; [2]->H2O; [1]->O2; [3]->OH
// Butler-Volmer prefactors: exchange current, active surface fraction and
// concentration-loss correction for the negative and positive electrodes.
//******************************************* constants used in eqation
	
		volScalarField C1Ne
		(
		Ne*as[0]*react.i0()[0]*(1-theta)*(CRa_Ne) 
		);

		volScalarField C1Pe
		(
		Pe*as[1]*react.i0()[1]*(1-theta)*(CRa_Pe) 
		);

		volScalarField C2Ne
		(
		Ne*as[0]*react.i0()[0]*(1-theta)*(CRc_Ne) 
		);

		volScalarField C2Pe
		(
		Pe*as[1]*react.i0()[1]*(1-theta)*(CRc_Pe) 
		);

// Exponential Butler-Volmer coefficients alpha*F/(R*T) for the anodic and
// cathodic branches at each electrode.
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

// First-order linearization of Butler-Volmer kinetics about the previous
// potential. These matrices provide implicit source terms for Ue/UNe/UPe.	
//*************************************** for the electrolyte
		
		fvScalarMatrix	Je_Ne
		(
		I1e_Ne*(1+C3Ne*Ue.oldTime())-fvm::Sp(I1e_Ne*C3Ne,Ue)
		+I2e_Ne*(-1-C4Ne*Ue.oldTime())+fvm::Sp(I2e_Ne*C4Ne,Ue)
		);
		fvScalarMatrix Je_Pe
		(
		I1e_Pe*(1+C3Pe*Ue.oldTime())-fvm::Sp(I1e_Pe*C3Pe,Ue)
		+I2e_Pe*(-1-C4Pe*Ue.oldTime())+fvm::Sp(I2e_Pe*C4Pe,Ue)
		);

//****************************************** for the electrodes
		
		fvScalarMatrix	J_Ne
		(
		I1_Ne*(1-C3Ne*UNe.oldTime())+fvm::Sp(I1_Ne*C3Ne,UNe)
		+I2_Ne*(-1+C4Ne*UNe.oldTime())+fvm::Sp(-I2_Ne*C4Ne,UNe)
		);
		fvScalarMatrix J_Pe
		(
		I1_Pe*(1-C3Pe*UPe.oldTime())+fvm::Sp(I1_Pe*C3Pe,UPe)
		+I2_Pe*(-1+C4Pe*UPe.oldTime())+fvm::Sp(-I2_Pe*C4Pe,UPe)
		);

        	
        // Effective electrolyte conductivity includes porous-medium and gas
        // saturation corrections consistent with paper Eq. 25.
        volScalarField sigmae
		(
		calKappa(T,C2[3]/1000)*(pow(epsilon1,poroM.tau())*pow(max(alpha2,SMALL),1.5))
		);
		
        	fvScalarMatrix UeEqn
        	(
        	fvm::laplacian(sigmae, Ue)
        	
        	==-1*(Je_Pe+Je_Ne)
        	);
		
		    // Solid-phase current continuity in the negative and positive
        	// porous electrodes (paper Eqs. 21-23).
        	//define potential equation for negative electrode.
        	
        	fvScalarMatrix UNeEqn
        	(fvm::laplacian(poroM.sigma_s_eff(), UNe)	==J_Ne	);
        	
        	//define potential equation for negative electrode.

        	fvScalarMatrix UPeEqn
        	(fvm::laplacian(poroM.sigma_s_eff(), UPe)	==J_Pe	);
        	
        	//solve eqations.
        	fvOptions.constrain(UNeEqn);
        	fvOptions.constrain(UeEqn);
        	fvOptions.constrain(UPeEqn);
            	UNeEqn.solve();

            	UeEqn.solve();

            	UPeEqn.solve();

            	
            	// Re-evaluate overpotential and volumetric transfer current with
            	// the newly solved potentials; J couples back to species production.
            	//evaluates J for the solving of C_[i] and plotting.
            	eta=(UNe-Ue-react.E0()[0])*Ne+(UPe-Ue-react.E0()[1])*Pe;
            	eta.correctBoundaryConditions();
        	J=
        	Pe*as[1]*react.i0()[1]*(1-theta)*(CRa_Pe*exp(C3Pe*eta)-(CRc_Pe)*exp(C4Pe*eta))+
        	Ne*as[0]*react.i0()[0]*(1-theta)*(CRa_Ne*exp(C3Ne*eta)-CRc_Ne*exp(C4Ne*eta))
    		;
		J.correctBoundaryConditions();
		Ie=fvc::reconstruct(UeEqn.flux());
