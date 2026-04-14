{
    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");
    Pair<tmp<volScalarField>> vDot_rhoC1 =
        mSTa.vDot_rhoC1();
    const volScalarField& vDotcAlphad = vDot_rhoC1[0]();
    const volScalarField& vDoteAlphad = vDot_rhoC1[1]();
    const volScalarField vDotvmcAlphad(vDotcAlphad - vDoteAlphad);
    
    tmp<surfaceScalarField> talphaPhi;
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
            fvm::Sp(vDotvmcAlphad, alpha1)
          + vDoteAlphad
        );

        alpha1Eqn.solve();

        Info<< "Phase-1 volume fraction = "
            << alpha1.weightedAverage(mesh.Vsc()).value()
            << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
            << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
            << endl;

        talphaPhi = alpha1Eqn.flux();
    }
    
    volScalarField alpha10("alpha10", alpha1);

    for (int aCorr=0; aCorr<nAlphaCorr; aCorr++)
    {
        tmp<surfaceScalarField> talphaPhiCorr
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
                alpha1,
                alpharScheme
            )
        );

        if (MULESCorr)
        {
            talphaPhiCorr.ref() -= talphaPhi();

            volScalarField alpha100("alpha100", alpha10);
            alpha10 = alpha1;

            MULES::correct
            (
                geometricOneField(),
                alpha1,
                talphaPhi(),
                talphaPhiCorr.ref(),
                (vDotvmcAlphad)(),
                (
                 - vDotvmcAlphad*alpha10
                )(),
                oneField(),
                zeroField()
            );

            // Under-relax the correction for all but the 1st corrector
            if (aCorr == 0)
            {
                talphaPhi.ref() += talphaPhiCorr();
            }
            else
            {
                alpha1 = 0.5*alpha1 + 0.5*alpha10;
                talphaPhi.ref() += 0.5*talphaPhiCorr();
            }
        }
        else
        {
            MULES::explicitSolve
            (
                geometricOneField(),
                alpha1,
                phi,
                talphaPhiCorr.ref(),
                (vDotvmcAlphad)(),
                (vDoteAlphad)(),
                oneField(),
                zeroField()
            );

            talphaPhi = talphaPhiCorr;
        }

        alpha2 = 1.0 - alpha1;
    }
    rhoPhi = talphaPhi()*(fvc::interpolate(mixture.rhod()) - mixture.rhoc()) + phi*mixture.rhoc();

    // Cache alphaPhi
    alphaPhi10 = talphaPhi();
    alphaPhi = talphaPhi();

    Info<< "Liquid phase volume fraction = "
        << alpha1.weightedAverage(mesh.V()).value()
        << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
        << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
        << endl;
}
