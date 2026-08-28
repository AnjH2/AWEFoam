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

#include "Dahlkild.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseDiffusionModels
{
    defineTypeNameAndDebug(Dahlkild, 0);
    addToRunTimeSelectionTable(phaseDiffusionModel, Dahlkild, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseDiffusionModels::Dahlkild::Dahlkild
(
    const dictionary& dict,
    const incompressibleTwoPhaseInteractingMixture& mixture,
    const word& modelName
)
:
    phaseDiffusionModel(dict, mixture,modelName),
    g_(meshObjects::gravity::New(mixture.U().time())),
    
    //dimensionedScalar
    Dp_("diffCoeffGravity", dimless, dict),
    Dn_("diffCoeffGravityNormal", dimless, dict),
    n_("n",dimless,dict),
    rb_("bubbleRadius",dimLength,dict),
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
        ),
    Solid_(
	    alphac_.mesh().lookupObject<volScalarField>
            (
            	"Solid"
            )
        )
{

    D_ =
        Dn_*tensor(1,0,0,0,1,0,0,0,1)
      + (Dp_ - Dn_)*(ep_*ep_);
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseDiffusionModels::Dahlkild::~Dahlkild()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
volScalarField Foam::phaseDiffusionModels::Dahlkild::kappa()
{
	
	return 0.6*pow(alphad_,2);
}

volScalarField Foam::phaseDiffusionModels::Dahlkild::beta()
{
	
	return (1.0/3.0)*pow(alphad_,2.0)*(1.0+0.5*exp(8.8*alphad_));
}

volScalarField Foam::phaseDiffusionModels::Dahlkild::f()
{
	
	return alphac_/(1/(alphac_+SMALL));
}

dimensionedVector Foam::phaseDiffusionModels::Dahlkild::vStokes()
{
	
	return g_.value()*(sqr(2*rb_)/(3*mixture_.nuc()))*dimensionedScalar(dimAcceleration,1);
}

volTensorField Foam::phaseDiffusionModels::Dahlkild::gamma()
{

    volTensorField gU_=fvc::grad(U_);
	return gU_;
}


//------	Diffusion coefficient functions 	------//

volTensorField Foam::phaseDiffusionModels::Dahlkild::UHdiff()
{
    // Terminal-velocity-based anisotropic dispersion contribution.

	return rb_*1/alphad_*f()*mag(vStokes())*D_;
}

volTensorField Foam::phaseDiffusionModels::Dahlkild::USdiff()
{

    // Isotropic velocity-gradient-driven dispersion contribution.
	//volTensorField	DOne_(dimless,tensor(1,0,0,0,1,0,0,0,1));	
	
	return tensor(1,0,0,0,1,0,0,0,1)*(pow(rb_,2)*mag(gamma())*beta()*1/alphad_);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::phaseDiffusionModels::Dahlkild::correct()
{
    // Convert the combined relative-diffusion closure to the mixture-mass
    // reference frame used by the solver.
    Ddm_=dF_*RToM()*(UHdiff()+USdiff());
 
}


// ************************************************************************* //
