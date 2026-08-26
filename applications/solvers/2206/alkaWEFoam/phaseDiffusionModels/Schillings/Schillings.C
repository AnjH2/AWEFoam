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
    DVp_(Dp_*ep_),
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

    	D_.replace(tensor::XX,mag(DVn1_.component(vector::X)+DVn2_.component(vector::X)+DVp_.component(vector::X)));
    	D_.replace(tensor::YY,mag(DVn1_.component(vector::Y)+DVn2_.component(vector::Y)+DVp_.component(vector::Y)));
    	D_.replace(tensor::ZZ,mag(DVn1_.component(vector::Z)+DVn2_.component(vector::Z)+DVp_.component(vector::Z)));
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
volScalarField Foam::phaseDiffusionModels::Schillings::kappa()
{
	
	return 0.6*pow(alphad_,2);
}

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
	
	return -g_.value()*(sqr(2.0*(rd_[0]*Ne_+rd_[1]*Pe_))/(18.0*mixture_.nuc()))*dimensionedScalar(dimAcceleration,1);
}

volScalarField Foam::phaseDiffusionModels::Schillings::gamma()
{

    //volTensorField tauP=dev2(fvc::grad(U_));
	//return fvc::grad(U_);//sqrt(0.5*( tauP && tauP ));
	//return sqrt(0.5*( tauP && tauP ));
	    volTensorField gU_=fvc::grad(U_);
   volTensorField E_=0.5*(gU_+gU_.T());
	return sqrt(2*( E_ && E_ ));
}


//------	Diffusion coefficient functions 	------//

volTensorField Foam::phaseDiffusionModels::Schillings::UHdiff()
{
    	//D_.replace(tensor::XX,DVn1_.component(vector::X)+DVn2_.component(vector::X)+DVp_.component(vector::X));
    	//D_.replace(tensor::YY,DVn1_.component(vector::Y)+DVn2_.component(vector::Y)+DVp_.component(vector::Y));
    	//D_.replace(tensor::ZZ,DVn1_.component(vector::Z)+DVn2_.component(vector::Z)+DVp_.component(vector::Z));	

	return 1/alphad_*(rd_[0]*Ne_+rd_[1]*Pe_)*f()*mag(vStokes())*D_;
}

volTensorField Foam::phaseDiffusionModels::Schillings::USdiff()
{
	//volTensorField	DOne_(dimless,tensor(1,0,0,0,1,0,0,0,1));
	return 1/alphad_*tensor(1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0)*pow((rd_[0]*(Ne_+NeC_)+rd_[1]*(Pe_+PeC_)),2)*mag(gamma())*beta();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::phaseDiffusionModels::Schillings::correct()
{
    Ddm_=dF_*RToM()*(UHdiff()+USdiff());

}


// ************************************************************************* //
