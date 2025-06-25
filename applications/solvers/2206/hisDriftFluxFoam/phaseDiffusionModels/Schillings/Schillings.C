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
    g_(meshObjects::gravity::New(mixture.U().time())),
    
    //dimensionedScalar
    Dp_("diffCoeffGravity", dimless, dict),
    Dn_("diffCoeffGravityNormal", dimless, dict),
    n_("n",dimless,dict),
    rb_("bubbleRadius",dimLength,dict),
ep_(vector(mag(g_.value().component(vector::X)),mag(g_.value().component(vector::Y)),mag(g_.value().component(vector::Z)))/mag(g_.value())),
    t_((mag(ep_.x()) < 0.9) ? vector(1,0,0) : vector(0,1,0)),
    en1_(ep_ ^ t_),
    en2_(ep_ ^ en1_),
    D_ 
    (
            IOobject
            (
                "anIsoDTensor",
                alphad_.mesh().time().timeName(),
                alphad_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE         // run-time field, not written
            ),
            alphad_.mesh(),
            dimensionedTensor("zero", dimless, tensor::zero) // start at 0
    ),
    DVn1_(Dn_*en1_),
    DVn2_(Dn_*en2_),
    DVp_(Dp_*ep_)
{

    	D_.replace(tensor::XX,mag(DVn1_.component(vector::X)+DVn2_.component(vector::X)+DVp_.component(vector::X)));
    	D_.replace(tensor::YY,mag(DVn1_.component(vector::Y)+DVn2_.component(vector::Y)+DVp_.component(vector::Y)));
    	D_.replace(tensor::ZZ,mag(DVn1_.component(vector::Z)+DVn2_.component(vector::Z)+DVp_.component(vector::Z)));
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseDiffusionModels::Schillings::~Schillings()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
volScalarField Foam::phaseDiffusionModels::Schillings::kappa()
{
	
	return 0.6*pow(alphad_,2);
}

volScalarField Foam::phaseDiffusionModels::Schillings::beta()
{
	
	return 1/3*pow(alphad_,2)*(1+0.5*exp(8.8*alphad_));
}

volScalarField Foam::phaseDiffusionModels::Schillings::f()
{
	
	return pow(1-alphad_,n_);
}

dimensionedVector Foam::phaseDiffusionModels::Schillings::vStokes()
{
	
	return -g_.value()*(sqr(2*rb_)/(18*mixture_.nuc()))*dimensionedScalar(dimAcceleration,1);
}

volScalarField Foam::phaseDiffusionModels::Schillings::gamma()
{

	return (en1_ & fvc::grad( mixture_.U() & ep_ ))+(en2_ & fvc::grad( mixture_.U() & ep_ ));
}


//------	Diffusion coefficient functions 	------//

volTensorField Foam::phaseDiffusionModels::Schillings::UHdiff()
{
    	//D_.replace(tensor::XX,DVn1_.component(vector::X)+DVn2_.component(vector::X)+DVp_.component(vector::X));
    	//D_.replace(tensor::YY,DVn1_.component(vector::Y)+DVn2_.component(vector::Y)+DVp_.component(vector::Y));
    	//D_.replace(tensor::ZZ,DVn1_.component(vector::Z)+DVn2_.component(vector::Z)+DVp_.component(vector::Z));	

	return rb_*f()*mag(vStokes())*D_;
}

volTensorField Foam::phaseDiffusionModels::Schillings::USdiff()
{
	//volTensorField	DOne_(dimless,tensor(1,0,0,0,1,0,0,0,1));
	return tensor(1,0,0,0,1,0,0,0,1)*pow(rb_,2)*(mag(gamma()))*beta();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::phaseDiffusionModels::Schillings::correct()
{
    Ddm_=(1-Mem_)*VToM()*(UHdiff()+USdiff());
}


// ************************************************************************* //
