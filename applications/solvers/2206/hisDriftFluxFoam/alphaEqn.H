{
    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");
    Pair<tmp<volScalarField>> vDotAlphal =
        mSTa.vDot_rhoC1();
    const volScalarField& vDotcAlphal = vDotAlphal[0]();
    const volScalarField& vDoteAlphal = vDotAlphal[1]();
    const volScalarField vDotvmcAlphal(vDotcAlphal - vDoteAlphal);
    
    if (MULESCorr)
    {
        fvScalarMatrix alpha1Eqn
        (
        fv::EulerDdtScheme<scalar>(mesh).fvmDdt(epsilon1,alpha1)
          + fv::gaussConvectionScheme<scalar>
            (
                mesh,
                phi,
                upwind<scalar>(mesh, phi)
            ).fvmDiv(phi, alpha1)
         ==
            fvm::Sp(vDotvmcAlphal, alpha1)
          + vDoteAlphal
        );

        solve(alpha1Eqn);

        Info<< "Phase-1 volume fraction = "
            << alpha1.weightedAverage(mesh.Vsc()).value()
            << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
            << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
            << endl;

        tmp<surfaceScalarField> talphaPhiUD(alpha1Eqn.flux());
        alphaPhi = talphaPhiUD();

        // Cache the upwind-flux
        talphaPhiCorr0 = talphaPhiUD;
    }

	volScalarField alpha10("alpha10", alpha1);
	
    for (int aCorr=0; aCorr<nAlphaCorr; aCorr++)
    {
    	
        tmp<surfaceScalarField> talphaPhiUn
        (
            fvc::flux
            (
                phi,
                alpha1,
                alphaScheme
            )
          + fvc::flux
            (
                phir,
                max(alpha1-dimensionedScalar(dimless,1e-4),Zero),
                alpharScheme
            )
        );

        if (MULESCorr)
        {
            tmp<surfaceScalarField> talphaPhiCorr(talphaPhiUn() - alphaPhi);
            alpha10 = alpha1;

            MULES::correct
            (
                epsilon1,
                alpha1,
                talphaPhiUn(),
                talphaPhiCorr.ref(),
                vDotvmcAlphal,
                (-vDotvmcAlphal*alpha1)(),
                UniformField<scalar>(mixture.alphaMax()),
                zeroField()
            );

            // Under-relax the correction for all but the 1st corrector
            if (aCorr == 0)
            {
                alphaPhi += talphaPhiCorr();
            }
            else
            {
                alpha1 = 0.5*alpha1 + 0.5*alpha10;
                alphaPhi += 0.5*talphaPhiCorr();
            }
        }
        else
        {
            alphaPhi = talphaPhiUn;
		
            MULES::explicitSolve
            (
                epsilon1,
                alpha1,
                phi,
                alphaPhi,
                vDotvmcAlphal,
                (vDoteAlphal)(),
                UniformField<scalar>(mixture.alphaMax()),
                zeroField()
            );
        }
    }

    alpha2 = 1.0 - alpha1;

    Info<< "Phase-1 volume fraction = "
        << alpha1.weightedAverage(mesh.Vsc()).value()
        << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
        << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
        << endl;
}
