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

#include "solidPressure.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseDiffusionModels
{
    defineTypeNameAndDebug(solidPressure, 0);
    addToRunTimeSelectionTable(phaseDiffusionModel, solidPressure, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseDiffusionModels::solidPressure::solidPressure
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
    A_ 
    (
            IOobject
            (
                "slopeGrowth",
                alphad_.mesh().time().timeName(),
                alphad_.mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE         // run-time field, not written
            ),
            alphad_.mesh(),
            dimensionedScalar("slopeGrowth", dimless, 1) // start at 0
    ),
    Cb_ 
    (
            IOobject
            (
                "diffusionScale",
                alphad_.mesh().time().timeName(),
                alphad_.mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE         // run-time field, not written
            ),
            alphad_.mesh(),
            dimensionedScalar("scale", dimless, 1) // start at 0
    ),
    alpha0_
    (
            IOobject
            (
                "diffusionOffset",
                alphad_.mesh().time().timeName(),
                alphad_.mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE         // run-time field, not written
            ),
            alphad_.mesh(),
            dimensionedScalar
            (
                "diffusionOffset",
                dimless,
                dict_.lookupOrDefault<scalar>("diffusionOffset", 0)
            )
    )
{

    	D_.replace(tensor::XX,mag(DVn1_.component(vector::X)+DVn2_.component(vector::X)+DVp_.component(vector::X)));
    	D_.replace(tensor::YY,mag(DVn1_.component(vector::Y)+DVn2_.component(vector::Y)+DVp_.component(vector::Y)));
    	D_.replace(tensor::ZZ,mag(DVn1_.component(vector::Z)+DVn2_.component(vector::Z)+DVp_.component(vector::Z)));
    if (!A_.headerOk())
    {
    A_ =
        dimensionedScalar("slopeGrowth_Ne",dimless,dict_)*(Ne_+NeC_)
      + dimensionedScalar("slopeGrowth_Pe",dimless,dict_)*(Pe_+PeC_)
      + dimensionedScalar("slopeGrowth_min", dimless, SMALL);
        
    }
    if (!Cb_.headerOk())
    {
    Cb_ =
        dimensionedScalar("diffusionScale_Ne",dimless,dict_)*(Ne_+NeC_+Mem_/2)
      + dimensionedScalar("diffusionScale_Pe",dimless,dict_)*(Pe_+PeC_+Mem_/2);
        
    }


}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseDiffusionModels::solidPressure::~solidPressure()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


volScalarField Foam::phaseDiffusionModels::solidPressure::f()
{
	return Cb_*exp(-1.0*A_*(alpha0_-alphad_));
}

volVectorField Foam::phaseDiffusionModels::solidPressure::vStokes()
{
	
	return -g_.value()*(sqr(2*rb_)/(18*mixture_.nuc()))*dimensionedScalar(dimAcceleration,1);
}



//------	Diffusion coefficient functions 	------//

volTensorField Foam::phaseDiffusionModels::solidPressure::UHdiff()
{
    	//D_.replace(tensor::XX,DVn1_.component(vector::X)+DVn2_.component(vector::X)+DVp_.component(vector::X));
    	//D_.replace(tensor::YY,DVn1_.component(vector::Y)+DVn2_.component(vector::Y)+DVp_.component(vector::Y));
    	//D_.replace(tensor::ZZ,DVn1_.component(vector::Z)+DVn2_.component(vector::Z)+DVp_.component(vector::Z));	

	return 1/alphad_*rb_*f()*mag(vStokes())*D_;
}




// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::phaseDiffusionModels::solidPressure::correct()
{
    Ddm_=dF_*(1-Mem_)*RToM()*UHdiff();
}


// ************************************************************************* //
