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

#include "BearScheidegger.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseDiffusionModels
{
    defineTypeNameAndDebug(BearScheidegger, 0);
    addToRunTimeSelectionTable(phaseDiffusionModel, BearScheidegger, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseDiffusionModels::BearScheidegger::BearScheidegger
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
    U_(
	        alphad_.mesh().lookupObject<volVectorField>
        	(
            		"U"
        	)
	),
    //dimensionedScalar
    Dp_("diffCoeffGravity", dimless, dict),
    Dn_("diffCoeffGravityNormal", dimless, dict),
    //n_("n",dimless,dict),
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
    aL_(dict.lookupOrDefault<scalar>("aL", 1)),
    bP_(dict.lookupOrDefault<scalar>("bP", 5)),
    cL_(dict.lookupOrDefault<scalar>("cL", 1)),
    dP_(dict.lookupOrDefault<scalar>("dP", 5)),
    
    BSPore0_("DispersitivCoefficient", dimLength, dict),
    nSat_(dict.lookupOrDefault<scalar>("nSat", -2)),
    chi_(dict.lookupOrDefault<scalar>("chi", 0.2))
    
{

    	D_.replace(tensor::XX,mag(DVn1_.component(vector::X)+DVn2_.component(vector::X)+DVp_.component(vector::X)));
    	D_.replace(tensor::YY,mag(DVn1_.component(vector::Y)+DVn2_.component(vector::Y)+DVp_.component(vector::Y)));
    	D_.replace(tensor::ZZ,mag(DVn1_.component(vector::Z)+DVn2_.component(vector::Z)+DVp_.component(vector::Z)));

}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseDiffusionModels::BearScheidegger::~BearScheidegger()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


volScalarField Foam::phaseDiffusionModels::BearScheidegger::f()
{
	
	return aL_*pow(1-alphad_,bP_)+cL_*pow(alphad_,dP_);
}

volVectorField Foam::phaseDiffusionModels::BearScheidegger::vStokes()
{
	
	return -g_.value()*(sqr(2*rb_)/(18*mixture_.nuc()))*dimensionedScalar(dimAcceleration,1);
}



//------	Diffusion coefficient functions 	------//

volTensorField Foam::phaseDiffusionModels::BearScheidegger::UHdiff()
{
    	//D_.replace(tensor::XX,DVn1_.component(vector::X)+DVn2_.component(vector::X)+DVp_.component(vector::X));
    	//D_.replace(tensor::YY,DVn1_.component(vector::Y)+DVn2_.component(vector::Y)+DVp_.component(vector::Y));
    	//D_.replace(tensor::ZZ,DVn1_.component(vector::Z)+DVn2_.component(vector::Z)+DVp_.component(vector::Z));	

	return 1/alphad_*rb_*f()*mag(vStokes())*D_;
}
volTensorField Foam::phaseDiffusionModels::BearScheidegger::UBSdiff()
{
    volVectorField  Uhat_=max(U_/eps_,dimensionedVector(dimVelocity,vector(SMALL,SMALL,SMALL)));
    volScalarField  BS_=BSPore0_*pow(alphac_,nSat_);

    return (
                        mag(Uhat_)*BS_*chi_*symmTensor::I
                        +
                        BS_*(1-chi_)*(Uhat_*Uhat_)/mag(Uhat_)
                    );
}



// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::phaseDiffusionModels::BearScheidegger::correct()
{

    Ddm_=dF_*(1-Mem_)*(RToM()*UHdiff()+alphac_/alphad_*VToM()*UBSdiff()/alphad_);

}


// ************************************************************************* //
