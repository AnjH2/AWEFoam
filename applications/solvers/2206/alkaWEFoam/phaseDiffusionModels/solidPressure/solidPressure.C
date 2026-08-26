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
    ANe_(dict.lookupOrDefault<scalar>("slopeGrowth_Ne", 1)),
    CbNe_(dict.lookupOrDefault<scalar>("scale_Ne", 1)),
    APe_(dict.lookupOrDefault<scalar>("slopeGrowth_Pe", 1)),
    CbPe_(dict.lookupOrDefault<scalar>("scale_Pe", 1)),
    alpha0_(dict.lookupOrDefault<scalar>("diffusionOffset", 0))
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

Foam::phaseDiffusionModels::solidPressure::~solidPressure()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


volScalarField Foam::phaseDiffusionModels::solidPressure::f()
{
	return CbNe_*exp(-1.0*ANe_*(alpha0_-alphad_))*(Ne_+NeC_+0.5*Mem_)+CbPe_*exp(-1.0*APe_*(alpha0_-alphad_))*(Pe_+PeC_+0.5*Mem_);
}

volVectorField Foam::phaseDiffusionModels::solidPressure::vStokes()
{
	
	return -g_.value()*(sqr(2*(rd_[0]*Ne_+rd_[1]*Pe_))/(18*mixture_.nuc()))*dimensionedScalar(dimAcceleration,1);
}



//------	Diffusion coefficient functions 	------//

volTensorField Foam::phaseDiffusionModels::solidPressure::UHdiff()
{
    	//D_.replace(tensor::XX,DVn1_.component(vector::X)+DVn2_.component(vector::X)+DVp_.component(vector::X));
    	//D_.replace(tensor::YY,DVn1_.component(vector::Y)+DVn2_.component(vector::Y)+DVp_.component(vector::Y));
    	//D_.replace(tensor::ZZ,DVn1_.component(vector::Z)+DVn2_.component(vector::Z)+DVp_.component(vector::Z));	

	return 1/alphad_*(rd_[0]*Ne_+rd_[1]*Pe_)*f()*mag(vStokes())*D_;
}




// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::phaseDiffusionModels::solidPressure::correct()
{
    Ddm_=dF_*(1-Mem_)*RToM()*UHdiff();
}


// ************************************************************************* //
