// Bounded gas-saturation transport. MULES advances the dispersed-phase
// fraction with mixture convection, relative/drift flux and interphase mass
// transfer while enforcing the prescribed saturation limits (paper Eq. 4).
{
    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");
    
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
            -mDote/rho1
        );

        solve(alpha1Eqn);

        Info<< "Phase-1 volume fraction = "
            << alpha1.weightedAverage(mesh.Vsc()).value()
            << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
            << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
            << endl;

        tmp<surfaceScalarField> talphaPhiUD(alpha1Eqn.flux());
        alphaPhi = talphaPhiUD();

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
                zeroField(),
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
                epsilon1,
                alpha1,
                phi,
                alphaPhi,
                zeroField(),
                (-mDote/rho1)(),
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
