{

//------------------Starting to correct permeability------------------

volVectorField Udm((1-Solid)*freeUdmModel.Udm()+Solid*UdmModel.Udm());
volVectorField UdmRd(Udm);
surfaceVectorField Udmf(fvc::interpolate(Udm));
/*const surfaceScalarField absDSolidf = Foam::mag(fvc::snGrad(Solid)*mesh.magSf());  // φ_nb - φ_P  with correct sign
// owner–neighbour addressing

volScalarField chi(fvc::surfaceSum(absDSolidf)*dimensionedScalar(dimless/dimLength,1));
// final normalisation χₚ = Σ_f |Δφ| / max(Σ_f |Δφ|, ξ)
forAll(chi, cellI)
{
    const scalar sum = chi[cellI];
    chi[cellI] = sum/(max(sum,VSMALL));          // VSMALL ≈ 1e-15
}

surfaceScalarField absOmegaf(Foam::mag(fvc::snGrad(chi*Solid)*mesh.magSf()));  // Δ(Δφ)_f  := dPhi_nb − dPhi_P
surfaceScalarField theta(absOmegaf);
forAll(absOmegaf, fI)
{
    theta[fI] = absOmegaf[fI]/(max(absOmegaf[fI], VSMALL));
}
// distance from cell centre to face centre (needed for ω_f)
surfaceScalarField w = 1.0/Foam::mag(mesh.delta());  // ω_f*/
//Info<< "sizes:  w1   = " << w1.size()   << nl   // should be 25275
//    << "        theta2 = " << theta2.size() << nl   // idem
//    << "        Udmf = " << Udmf.size() << endl;    // idem
//Info<<"did we get this far"<<endl;
forAll(activeCells2, idx)
{
    //if (chi2[cellI] < 0.5) continue;     // keep original value
    const label cellI = activeCells2[idx];
    vector num = vector::zero;          // numerator  (vector)
    scalar den = 0.0;                   // denominator (scalar)
    //Info<<"idx:"<<idx<<"    cellI:"<<cellI<<endl;
    const labelUList& cFaces = mesh.cells()[cellI];

    
        // walk over faces of cellI
    forAll(cFaces, i)
    {
        //Info<<"cFaces:"<<cFaces[i]<<"    i:"<<i<<"  chi2[cellI]:"<<chi2[cellI]<<endl;
         const label fI = cFaces[i];
        if (fI >= mesh.nInternalFaces()) continue;
        scalar wgt = w1[fI] * theta2[fI];
        num += wgt * Udmf[fI];
        den += wgt;
    }
     //Info<<"num:"<<num<<"    den:"<<den<<endl;
    if (den > VSMALL)
        UdmRd[cellI] = num / den;
}
//Info<<"did we complete??"<<endl;
UdmRd.correctBoundaryConditions();
phir=fvc::flux(UdmRd);
}
