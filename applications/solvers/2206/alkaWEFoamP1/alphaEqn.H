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
            epsilon1*fv::EulerDdtScheme<scalar>(mesh).fvmDdt(alpha1)
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
                phi*fvc::interpolate(1/epsilon1),
                alpha1,
                alphaScheme
            )
          + fvc::flux
            (
                phir,
                alpha1,
                alpharScheme
            )
        );

        if (MULESCorr)
        {
            tmp<surfaceScalarField> talphaPhiCorr(talphaPhiUn() - alphaPhi);
            volScalarField alpha100("alpha100", alpha10);
            alpha10 = alpha1;

            MULES::correct
            (
                geometricOneField(),
                alpha1,
                talphaPhiUn(),
                talphaPhiCorr.ref(),
                vDotvmcAlphal,
                (
                    (alpha10 - alpha100)
                  - vDotvmcAlphal*alpha10
                )(),
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
                geometricOneField(),
                alpha1,
                phi*fvc::interpolate(1/epsilon1),
                alphaPhi,
                vDotvmcAlphal,
                vDoteAlphal,
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
