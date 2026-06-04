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

#include "darcyForchheimer.H"
#include "addToRunTimeSelectionTable.H"
#include "../../porousProperties/porousProperties.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace relativeVelocityModels
{
    defineTypeNameAndDebug(darcyForchheimer, 0);
    addToRunTimeSelectionTable(relativeVelocityModel, darcyForchheimer, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::relativeVelocityModels::darcyForchheimer::darcyForchheimer
(
    const dictionary& dict,
    const incompressibleTwoPhaseInteractingMixture& mixture,
    const word& modelName
)
:
    relativeVelocityModel(dict, mixture,modelName),
    mixture_(mixture),
    electrodes({"Ne","Pe"}),
    dict_(dict),
    g_(meshObjects::gravity::New(mixture.U().time())),
    
    rb_
    (
        mixture_.alpha1().mesh().lookupObject<volScalarField>
        (
            "rb"
        )
    ),
    CF_(electrodes.size()),
	//(
	//	"R_DB",
	//	dimensionSet ( 0, 1, 0, 0, 0, 0,0),
	//	dict
	//),
    n_("n",dimless,dict_),
    minAlphad_("minAlphad",dimless,dict_),
    rhoc_(mixture.rhoc()),
    rhod_(mixture.rhod()),
    
    
    eps_(
	    alphac_.mesh().lookupObject<volScalarField>
            (
            	"eps"
            )
        ),
    p_(
	    alphac_.mesh().lookupObject<volScalarField>
            (
            	"p"
            )
        ),
    K_(
        eps_.mesh().lookupObject<porousProperties>
        (
            "porousProperties"
        ).K()
        )
{

}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::relativeVelocityModels::darcyForchheimer::~darcyForchheimer()
{}

// * * * * * * * * * * * * * * Private Functions  * * * * * * * * * * * * * * //
volScalarField Foam::relativeVelocityModels::darcyForchheimer::f()
{
	
	return pow(1-min(alphad_,minAlphad_),n_);
}
volScalarField Foam::relativeVelocityModels::darcyForchheimer::krd()
{
	
	return pow(alphad_,mixture_.kr());
}
volScalarField Foam::relativeVelocityModels::darcyForchheimer::krc()
{
	
	return pow(1-alphad_,mixture_.kr());
}

volScalarField Foam::relativeVelocityModels::darcyForchheimer::A()
{

	return 1/K_*(mixture_.mud_m()/(krd()+SMALL)+mixture_.muc()/(krc()+SMALL));
}



volScalarField Foam::relativeVelocityModels::darcyForchheimer::B()
{

	return ((PeC_+Pe_+Mem_*1/2)*CF_[1]+(NeC_+Ne_+Mem_*1/2)*CF_[0])*(rhoc_/sqrt(K_*(krc()+SMALL))+rhod_/sqrt(K_*(krd()+SMALL)));
}
volVectorField Foam::relativeVelocityModels::darcyForchheimer::force()
{
    //surfaceScalarField dpcdsf(fvc::interpolate(mixture_.dpcds(), "harmonic"));
    
	//return (rhod_-rhoc_)*g_-fvc::reconstruct(fvc::snGrad(mixture_.pc())* alphad_.mesh().magSf());
	
	return (rhoc_-rhod_)*g_-fvc::grad(p_);//-fvc::grad(mixture_.pc());
	
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::relativeVelocityModels::darcyForchheimer::correct()
{    

    dModel_->correct();


    tmp<volScalarField> u_=(-A()+sqrt(sqr(A())+4*B()*mag(force())))/(2*B());


    Udm_ = ((-1*rhoc_*(1-alphad_))/rho())*u_*force()/mag(force());
}


// ************************************************************************* //
