/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2015 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "Schillings.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseDiffusionModels
{
    defineTypeNameAndDebug(Schillings, 0);
    addToRunTimeSelectionTable(phaseDiffusionModel, Schillings, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseDiffusionModels::Schillings::Schillings
(
    const dictionary& dict,
    const incompressibleTwoPhaseInteractingMixture& mixture,
    const word& modelName
)
:
    phaseDiffusionModel(dict, mixture,modelName),
    electrodes({"Ne","Pe"}),
    dict_(dict),
    g_(meshObjects::gravity::New(mixture.U().time())),
    
    //dimensionedScalar
    Dp_("diffCoeffGravity", dimless, dict),
    Dn_("diffCoeffGravityNormal", dimless, dict),
    n_("n",dimless,dict),
    rd_(electrodes.size()),
    ep_(g_.value()/mag(g_.value())),

    D_
    (
        IOobject
        (
            "anIsoDTensor",
            alphad_.mesh().time().timeName(),
            alphad_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alphad_.mesh(),
        dimensionedTensor("zero", dimless, tensor::zero)
    ),
    U_(
	    alphac_.mesh().lookupObject<volVectorField>
            (
            	"U"
            )
        )
{
    D_ =
        Dn_*tensor(1,0,0,0,1,0,0,0,1)
      + (Dp_ - Dn_)*(ep_*ep_);
    forAll(electrodes,i)
	{
	rd_.set
    	(
        	i,
        	new dimensionedScalar("bubbleRadius_"+electrodes[i], dimLength,dict_)
        );
    }
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseDiffusionModels::Schillings::~Schillings()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


volScalarField Foam::phaseDiffusionModels::Schillings::beta()
{
	
	return (1.0/3.0)*pow(alphad_,2.0)*(1.0+0.5*exp(8.8*alphad_));
}

volScalarField Foam::phaseDiffusionModels::Schillings::f()
{
	
	return pow(1-alphad_,n_);
}

volVectorField Foam::phaseDiffusionModels::Schillings::vStokes()
{
	
	return -g_.value()*sqr(2.0*(rd_[0]*(Ne_+NeC_)+rd_[1]*(Pe_+PeC_)))/(18.0*mixture_.nuc())*dimensionedScalar(dimAcceleration,1);

}

volScalarField Foam::phaseDiffusionModels::Schillings::gamma()
{

	    volTensorField gU_=fvc::grad(U_);
   volTensorField E_=0.5*(gU_+gU_.T());
	return sqrt(2*( E_ && E_ ));
}


//------	Diffusion coefficient functions 	------//

volTensorField Foam::phaseDiffusionModels::Schillings::UHdiff()
{
	
    // Hydrodynamic self-diffusion scaled with bubble radius, terminal
    // velocity and the prescribed anisotropy tensor.
	return 1/alphad_*(rd_[0]*(Ne_+NeC_)+rd_[1]*(Pe_+PeC_))*f()*mag(vStokes())*D_;
}

volTensorField Foam::phaseDiffusionModels::Schillings::USdiff()
{
    // Isotropic shear-induced dispersion proportional to r_b^2 |gamma|
    // and the gas-fraction function beta().
	//volTensorField	DOne_(dimless,tensor(1,0,0,0,1,0,0,0,1));
	return 1/alphad_*tensor(1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0)*pow((rd_[0]*(Ne_+NeC_)+rd_[1]*(Pe_+PeC_)),2)*mag(gamma())*beta();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::phaseDiffusionModels::Schillings::correct()
{
    // Combine hydrodynamic and shear-induced diffusion independently of
    // the buoyancy/lift drift closure.
    Ddm_=dF_*RToM()*(UHdiff()+USdiff());

}


// ************************************************************************* //
