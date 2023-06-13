{
    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");

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
          
        );
	alpha1Eqn += Psi_gas;
	
        solve(alpha1Eqn);

        Info<< "Phase-1 volume fraction = "
            << alpha1.weightedAverage(mesh.Vsc()).value()
            << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
            << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
            << endl;

        tmp<surfaceScalarField> talphaPhiUD(alpha1Eqn.flux());
        alphaPhi = talphaPhiUD();

        if (alphaApplyPrevCorr && talphaPhiCorr0.valid())
        {
            Info<< "Applying the previous iteration correction flux" << endl;

            MULES::correct
            (
                geometricOneField(),
                alpha1,
                alphaPhi,
                talphaPhiCorr0.ref(),
                Psi_gas,
                zeroField(),
                UniformField<scalar>(mixture.alphaMax()),
                zeroField()
            );

            alphaPhi += talphaPhiCorr0();
        }

        // Cache the upwind-flux
        talphaPhiCorr0 = talphaPhiUD;
    }

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
                alpha1,
                alpharScheme
            )
        );

        if (MULESCorr)
        {
            tmp<surfaceScalarField> talphaPhiCorr(talphaPhiUn() - alphaPhi);
            volScalarField alpha10("alpha10", alpha1);

            MULES::correct
            (
                geometricOneField(),
                alpha1,
                talphaPhiUn(),
                talphaPhiCorr.ref(),
                Psi_gas,
                zeroField(),
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
                phi,
                alphaPhi,
                Psi_gas,
                zeroField(),
                UniformField<scalar>(mixture.alphaMax()),
                zeroField()
            );
        }
    }

    if (alphaApplyPrevCorr && MULESCorr)
    {
        talphaPhiCorr0 = alphaPhi - talphaPhiCorr0;
    }

    alpha2 = 1.0 - alpha1;

    Info<< "Phase-1 volume fraction = "
        << alpha1.weightedAverage(mesh.Vsc()).value()
        << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
        << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
        << endl;
        
    //try calculating alphaphi2 (continous phase)

    alpha2Phi=
        (
            fvc::flux
            (
                phi,
                alpha2,
                alphaScheme
            )
          + fvc::flux
            (
                -fvc::flux(phir, alpha1, alpharScheme),
                alpha2,
                alpharScheme
            )
        );
}
