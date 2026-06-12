
    phir=fvc::flux((1-Solid)*freeUdmModel.Udm()+(Solid)*UdmModel.Udm());
    //#include "correctUdm.H" 
     
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

    if (UdmModel.runDiffusion() or freeUdmModel.runDiffusion()) {
        Info<<"Running phase diffusion equation"<<endl;
        volSymmTensorField Ddm("Ddm",symm(Solid*UdmModel.Ddm()+(1-Solid)*freeUdmModel.Ddm()));
        surfaceScalarField nDn
        (
            "nDn",
            (mesh.Sf() & (fvc::interpolate(Ddm) & mesh.Sf()))
            / (sqr(mesh.magSf()) + VSMALL*dimensionedScalar(dimArea*dimArea,1))
        );
        
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

            
    }
    
    rhoPhi = alphaPhi*fvc::interpolate(mixture.rhod() - mixture.rhoc()) + phi*mixture.rhoc();
    rho = mixture.rho()*poroM.eps();
    rho.correctBoundaryConditions();
        surfaceScalarField alpha2Phi
        (
            //(rhoPhi-alphaPhi*fvc::interpolate(mixture.rhod()))/mixture.rhoc()
            (phi*fvc::interpolate(rho) - fvc::interpolate(epsilon1*mixture.rhod())*alphaPhi)/mixture.rhoc()
        );
    Ur=mixture.rho()/(alpha2*mixture.rhoc())*((1-Solid)*freeUdmModel.Udm()+Solid*UdmModel.Udm()-((Solid*epsilon1*UdmModel.Ddm()+(1-Solid)*freeUdmModel.Ddm()) & fvc::grad(alpha1)));

