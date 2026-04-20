{

//------------------Starting to correct permeability------------------

volScalarField D(mixture.mu()/poroM.K());
Drd=D;
surfaceScalarField Df(fvc::interpolate(D));
/*const surfaceScalarField absDEpsf = Foam::mag(fvc::snGrad(poroM.eps())*mesh.magSf());  // φ_nb - φ_P  with correct sign
// owner–neighbour addressing

volScalarField chi(fvc::surfaceSum(absDEpsf)*dimensionedScalar(dimless/dimLength,1));
// final normalisation χₚ = Σ_f |Δφ| / max(Σ_f |Δφ|, ξ)
forAll(chi, cellI)
{
    const scalar sum = chi[cellI];
    chi[cellI] = sum/(max(sum,VSMALL));          // VSMALL ≈ 1e-15
}

surfaceScalarField absOmegaf(Foam::mag(fvc::snGrad(chi*poroM.eps())*mesh.magSf()));  // Δ(Δφ)_f  := dPhi_nb − dPhi_P
surfaceScalarField theta(absOmegaf);
forAll(absOmegaf, fI)
{
    theta[fI] = absOmegaf[fI]/(max(absOmegaf[fI], VSMALL));
}
// distance from cell centre to face centre (needed for ω_f)
surfaceScalarField w = 1.0/Foam::mag(mesh.delta());  // ω_f

forAll(Drd, cellI)
{
    if (chi1[cellI] < 0.5) continue;                  // keep original D_P

    scalar num = 0.0, den = 0.0;

    // walk over faces of cellI
    forAll(mesh.Cf(), fI)
    {
        if (mesh.faceOwner()[fI] != cellI && mesh.faceNeighbour()[fI] != cellI) continue;

        num += w1[fI]*theta1[fI]*Df[fI];
        den += w1[fI]*theta1[fI];
    }

    Drd[cellI] = (den > VSMALL) ? num/den : D[cellI];
}
forAll(activeCells1, idx)
{
    //if (chi2[cellI] < 0.5) continue;     // keep original value
    const label cellI = activeCells1[idx];
    scalar num = 0.0, den = 0.0;
    //Info<<"idx:"<<idx<<"    cellI:"<<cellI<<endl;
    const labelUList& cFaces = mesh.cells()[cellI];

    
        // walk over faces of cellI
    forAll(cFaces, i)
    {
        //Info<<"cFaces:"<<cFaces[i]<<"    i:"<<i<<"  chi2[cellI]:"<<chi2[cellI]<<endl;
         const label fI = cFaces[i];
        if (fI >= mesh.nInternalFaces()) continue;
        scalar wgt = w1[fI] * theta1[fI];
        num += wgt * Df[fI];
        den += wgt;
    }
     //Info<<"num:"<<num<<"    den:"<<den<<endl;
    if (den > VSMALL)
        Drd[cellI] = (den > VSMALL) ? num/den : D[cellI];
}*/

volScalarField numP (fvc::surfaceSum(wgtf * Df));
forAll(activeCells1, idx)
{
    const label cellI = activeCells1[idx];
    Drd[cellI] = numP[cellI] / denP[cellI];
}



Drd.correctBoundaryConditions();

}
