
    //phir=fvc::flux((1-Solid)*freeUdmModel.Udm()+Solid*UdmModel.Udm());
    #include "correctUdm.H" 
      Pair<tmp<volScalarField>> mDot =
        mSTa.mDot();
    const volScalarField& mDotc = mDot[0]();
    const volScalarField& mDote = mDot[1]();

    surfaceScalarField alphaPhi
    (
        IOobject
        (
            "alphaPhi",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar(phi.dimensions(), Zero)
    );

    if (nAlphaSubCycles > 1)
    {
        dimensionedScalar totalDeltaT = runTime.deltaT();
        surfaceScalarField alphaPhiSum
        (
            IOobject
            (
                "alphaPhiSum",
                runTime.timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar(phi.dimensions(), Zero)
        );

        for
        (
            subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
            !(++alphaSubCycle).end();
        )
        {
            #include "alphaEqn.H"
            alphaPhiSum += (runTime.deltaT()/totalDeltaT)*alphaPhi;
        }

        alphaPhi = alphaPhiSum;
    }
    else
    {
        #include "alphaEqn.H"
    }
    

  //for (int aCorr=0; aCorr<nAlphaCorr; aCorr++)
  
	//loop to update mass transfer
	//for(label i=0; i<mixture.species1L().size(); i++)
	//{
	//	mSTa.correct_mDot_wall(i,sTp.C2_s(),theta);
        //}
        /*fvScalarMatrix alpha1Eqn
        (
            fvm::ddt(epsilon1*rho1,alpha1)
            //+fvm::div(phi*fvc::interpolate((1-Mem+VSMALL)*rho1),alpha1,"div(U,alphad)")
            +fvm::div(phi*fvc::interpolate((1-Mem+VSMALL)*rho1),alpha1,"div(U,alphad)")
            ==
           // -(
           // 	mSTa.mDotAlphal()[1]*pow(alpha1,2)+fvm::Sp(mSTa.mDotAlphal()[1]*(-2*alpha1+1), alpha1)
           // )//+mDotc
            -mSTa.mDot()[1]
            +fvm::laplacian((1-Mem+VSMALL)*rho1*UdmModel.Ddm(), alpha1)
            -fvm::div(phir*fvc::interpolate((1-Mem+VSMALL)*rho1),alpha1,"div(Udm,alphad)")
        );
        fvScalarMatrix alpha1Eqn
        (
            fvm::ddt(epsilon1*rho1,alpha1)
            //+fvm::div(phi*fvc::interpolate((1-Mem+VSMALL)*rho1),alpha1,"div(U,alphad)")
            +fvm::div(phi*fvc::interpolate(rho1),alpha1,"div(U,alphad)")
            ==
           // -(
           // 	mSTa.mDotAlphal()[1]*pow(alpha1,2)+fvm::Sp(mSTa.mDotAlphal()[1]*(-2*alpha1+1), alpha1)
           // )//+mDotc
            -mSTa.mDot()[1]
            +fvm::laplacian((1-Mem)*epsilon1*rho1*UdmModel.Ddm(), alpha1)
            -fvm::div(phir*fvc::interpolate((1-Mem+VSMALL)*rho1),alpha1,"div(Udm,alphad)")
        );
        alpha1Eqn.solve(mesh.solver("alpha1Eqn"));
	
        alphaPhi = alpha1Eqn.flux()/fvc::interpolate(rho1);
        
        alpha2 = 1.0 - alpha1;
        
        
	
	
        Info<< "Phase-1 volume fraction = "
            << alpha1.weightedAverage(mesh.Vsc()).value()
            << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
            << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
            << endl;
    }
	//alphaPhi=alphaPhi*fvc::interpolate(epsilon1); //to make it superficial
    //rhoPhi = alphaPhi*(fvc::interpolate(mixture.rhod()) - rho2) + phi*rho2; //original
    //rhoPhi = fvc::interpolate(alpha1*mixture.rhod())*(phi+phir)+fvc::interpolate(1-alpha1)*(phi+phiL)*rho2; //V665
    //rhoPhi = alphaPhi*(fvc::interpolate(mixture.rhod()) - rho2)+fvc::interpolate(1-alpha1)*(phi+phiL)*rho2; //V722
    //surfaceScalarField alpha2Phi=fvc::interpolate(1-alpha1)*(fvc::interpolate(1/alpha1)*alphaPhi-phir-phiDdm+phiDcm+phiL);
    //rhoPhi = fvc::interpolate(rho1)*alphaPhi+fvc::interpolate(1-alpha1)*(fvc::interpolate(1/alpha1)*alphaPhi-phir-phiDdm+phiDcm+phiL)*rho2; //V754
    rhoPhi = alphaPhi*fvc::interpolate(mixture.rhod() - mixture.rhoc()) + phi*fvc::interpolate(mixture.rhod());
    rho = mixture.rho();
    */
    //rhoPhi = fvc::interpolate(rho)*(fvc::interpolate(1/alpha1)*alphaPhi-phir);; //V755
    if (UdmModel.runDiffusion() or freeUdmModel.runDiffusion()) {
        Info<<"Running phase diffusion equation"<<endl;
        fvScalarMatrix alpha1Eqn
        (
            fvm::ddt(alpha1) - fvc::ddt(alpha1)
          - fvm::laplacian(turbulence->nut(), alpha1)
          - fvm::laplacian((1-Mem)*epsilon1*(Solid*UdmModel.Ddm()+(1-Solid)*freeUdmModel.Ddm()), alpha1)
        );

        alpha1Eqn.solve(mesh.solver("alpha1Diffusion"));
        alphaPhi += alpha1Eqn.flux();
        alpha2 = 1.0 - alpha1;

        Info<< "Phase-1 volume fraction = "
            << alpha1.weightedAverage(mesh.Vsc()).value()
            << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
            << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
            << endl;
    }

    rhoPhi = alphaPhi*fvc::interpolate(mixture.rhod() - mixture.rhoc()) + phi*fvc::interpolate(mixture.rhod());
    rho = mixture.rho()*poroM.eps();

