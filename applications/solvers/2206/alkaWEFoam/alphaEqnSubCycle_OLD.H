
    phir=fvc::flux((1-Solid)*freeUdmModel.Udm()+Solid*UdmModel.Udm());
    //#include "correctUdm.H" 
    divU=fvc::div(phi);
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
 /*
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
        Info<<"before alphaEqn"<<endl;
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
    Info<<"after alphaEqn"<<endl;
*/
      Pair<tmp<volScalarField>> mDot =
        mSTa.mDot();
    const volScalarField& mDotc = mDot[0]();
    const volScalarField& mDote = mDot[1]();
    //Info<<"mass Transfer courant number >> "<< max(mag(mDot[1]()*runTime.deltaTValue()/(epsilon1*rho1))).value()<<endl;
  for (int aCorr=0; aCorr<nAlphaCorr; aCorr++)
  {

	//loop to update mass transfer
	//for(label i=0; i<mixture.species1L().size(); i++)
	//{
	//	mSTa.correct_mDot_wall(i,sTp.C2_s(),theta);
        //}
        volSymmTensorField Ddm("Ddm",symm(Solid*UdmModel.Ddm()+(1-Solid)*freeUdmModel.Ddm()));
        surfaceScalarField nDn
        (
            "nDn",
            (mesh.Sf() & (fvc::interpolate(Ddm) & mesh.Sf()))
            / (sqr(mesh.magSf()) + VSMALL*dimensionedScalar(dimArea*dimArea,1))
        );
        //Info<<"---CHECK---      min(Ddm)->"<<min(Ddm)<<endl;
        fvScalarMatrix alpha1Eqn
        (
            fvm::ddt(epsilon1,alpha1)
            //+fvm::div(phi*fvc::interpolate((1-Mem+VSMALL)*rho1),alpha1,"div(U,alphad)")
            +fvm::div(phi,alpha1,"div(U,alphad)")
            
            ==
           // -(
           // 	mSTa.mDotAlphal()[1]*pow(alpha1,2)+fvm::Sp(mSTa.mDotAlphal()[1]*(-2*alpha1+1), alpha1)
           // )//+mDotc
            -mSTa.mDot()[1]/rho1
            -fvc::div(phir,max(alpha1-dimensionedScalar(dimless,1e-4),Zero),"div(Udm,alphad)")//+fvc::div(phir,max(alpha1-dimensionedScalar(dimless,1e-4),Zero),"div(Udm,alphad)")
            + fvm::laplacian(nDn*fvc::interpolate(epsilon1), alpha1,"laplacian(D,alpha1)")
            //- fvc::laplacian(nDn*fvc::interpolate(epsilon1),1e-4,"laplacian(D,alpha1)")
           
            
        );
        alpha1Eqn.solve(mesh.solver("alpha1Eqn"));
	
        alphaPhi = alpha1Eqn.flux();
        
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
    //rhoPhi = alphaPhi*fvc::interpolate(mixture.rhod() - mixture.rhoc()) + phi*fvc::interpolate(mixture.rhod());
    //rho = mixture.rho();
    
    //rhoPhi = fvc::interpolate(rho)*(fvc::interpolate(1/alpha1)*alphaPhi-phir);; //V755
    
    
    //2025-12-04 works well
    /*if (UdmModel.runDiffusion() or freeUdmModel.runDiffusion()) {
        Info<<"Running phase diffusion equation"<<endl;
        volSymmTensorField Ddm("Ddm",symm(Solid*UdmModel.Ddm()+(1-Solid)*freeUdmModel.Ddm()));
        surfaceScalarField nDn
        (
            "nDn",
            (mesh.Sf() & (fvc::interpolate(Ddm) & mesh.Sf()))
            / (sqr(mesh.magSf()) + VSMALL*dimensionedScalar(dimArea*dimArea,1))
        );
        Info<<"---CHECK---      min(Ddm)->"<<min(Ddm)<<endl;
        fvScalarMatrix alpha1Eqn
        (
            fvm::ddt(epsilon1,alpha1) - fvc::ddt(epsilon1,alpha1)
          //- fvm::laplacian(turbulence->nut(), alpha1)
          - fvm::laplacian(nDn*fvc::interpolate(epsilon1), alpha1,"laplacian(D,alpha1)")
        );

        alpha1Eqn.solve(mesh.solver("alpha1Diffusion"));
        alphaPhi += alpha1Eqn.flux();
        alpha2 = 1.0 - alpha1;

        Info<< "Phase-1 volume fraction = "
            << alpha1.weightedAverage(mesh.Vsc()).value()
            << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
            << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
            << endl;

            
    }*/
    
    rhoPhi = alphaPhi*fvc::interpolate(mixture.rhod() - mixture.rhoc()) + phi*fvc::interpolate(mixture.rhod());
    rho = mixture.rho()*poroM.eps();
    rho.correctBoundaryConditions();
        surfaceScalarField alpha2Phi
        (
            //(rhoPhi-alphaPhi*fvc::interpolate(mixture.rhod()))/mixture.rhoc()
            (phi*fvc::interpolate(rho) - fvc::interpolate(epsilon1*mixture.rhod())*alphaPhi)/mixture.rhoc()
        );
    Ur=mixture.rho()/(alpha2*mixture.rhoc())*((1-Solid)*freeUdmModel.Udm()+Solid*UdmModel.Udm()-((Solid*epsilon1*UdmModel.Ddm()+(1-Solid)*freeUdmModel.Ddm()) & fvc::grad(alpha1)));

