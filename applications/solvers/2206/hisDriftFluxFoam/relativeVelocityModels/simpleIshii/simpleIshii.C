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

#include "simpleIshii.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace relativeVelocityModels
{
    defineTypeNameAndDebug(simpleIshii, 0);
    addToRunTimeSelectionTable(relativeVelocityModel, simpleIshii, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::relativeVelocityModels::simpleIshii::simpleIshii
(
    const dictionary& dict,
    const incompressibleTwoPhaseInteractingMixture& mixture
)
:
    relativeVelocityModel(dict, mixture),
    mixture_(mixture),
    g_(meshObjects::gravity::New(mixture.U().time())),
    
    rd_
	(
		"R_DB",
		dimensionSet ( 0, 1, 0, 0, 0, 0,0),
		dict
	),
    n_("n",dimless,dict),
    rhoc_(mixture.rhoc()),
    rhod_(mixture.rhod()),
    D_("D", dimVelocity*dimLength, dict),
        V1_(
	alphac_.mesh().lookupObject<volScalarField>
        (
            	"V1"
        )
        ),
    eps_(
	    alphac_.mesh().lookupObject<volScalarField>
            (
            	"eps"
            )
        ),
    U_(
	    alphac_.mesh().lookupObject<volVectorField>
            (
            	"U"
            )
        ),
    eg_("eg",(-1*g_)/mag(g_)),
    en_(dict.lookupOrDefault<dimensionedVector>("en", dimensionedVector("en", dimless, vector (0,1,0)))),
    em_(dict.lookupOrDefault<dimensionedVector>("em", dimensionedVector("em", dimless, vector (0,0,1))))
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::relativeVelocityModels::simpleIshii::~simpleIshii()
{}

// * * * * * * * * * * * * * * Private Functions  * * * * * * * * * * * * * * //
volScalarField Foam::relativeVelocityModels::simpleIshii::kappa()
{
	
	return 0.6*pow(alphad_,2);
}

volScalarField Foam::relativeVelocityModels::simpleIshii::beta()
{
	
	return 1/3*pow(alphad_,2)*(1+0.5*exp(8.8*alphad_));
}

volScalarField Foam::relativeVelocityModels::simpleIshii::f()
{
	
	return pow(1-alphad_,n_);
}

volVectorField Foam::relativeVelocityModels::simpleIshii::vStokes()
{
	
	return mag(g_)*sqr(2*rd_)/(18*mixture_.nucModel().nu())*eg_;
}

volScalarField Foam::relativeVelocityModels::simpleIshii::gamma(const dimensionedVector i)
{
	
	
	return fvc::grad(U_&eg_)&i;
}

volVectorField Foam::relativeVelocityModels::simpleIshii::UStokes()
{
	
	return f()*mag(vStokes())*eg_*(1-Mem_+VSMALL);
}

volVectorField Foam::relativeVelocityModels::simpleIshii::USaff()
{
	
	return f()*mag(vStokes())*sign(gamma(en_))*6.46/(6*Foam::constant::mathematical::pi)*sqrt((sqr(rd_)*mag(gamma(en_)))/mixture_.nucModel().nu())*en_*(1-Mem_+VSMALL)
	 +f()*mag(vStokes())*sign(gamma(em_))*6.46/(6*Foam::constant::mathematical::pi)*sqrt((sqr(rd_)*mag(gamma(em_)))/mixture_.nucModel().nu())*em_*(1-Mem_+VSMALL);
}

volVectorField Foam::relativeVelocityModels::simpleIshii::USmig()
{
	return -1/alphad_*pow(rd_,2)*mag(gamma(en_))*kappa()/(mixture_.mu()*(dimensionedScalar(dimless/dimTime,VSMALL)+(fvc::grad(U_)&eg_&en_)))*fvc::grad(fvc::grad(U_*mixture_.mu())&eg_&en_)*(1-Mem_+VSMALL)
	-1/alphad_*pow(rd_,2)*mag(gamma(em_))*kappa()/(mixture_.mu()*(dimensionedScalar(dimless/dimTime,VSMALL)+(fvc::grad(U_)&eg_&em_)))*fvc::grad(fvc::grad(U_*mixture_.mu())&eg_&em_)*(1-Mem_+VSMALL);
}
// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::relativeVelocityModels::simpleIshii::correct()
{
    
    
    Udm_ = (1-Mem_+VSMALL)*(rhoc_/rho())*(UStokes()+USaff()+USmig());
    /*
    Udm_.component(0) = pow(2,0.5)*pow((sigma_*g_.component(0)*(rhod_-rhoc_))/(pow(rhoc_,2)),0.25)*pow(1-alphad_,1.75);
    Info<<Udm_.component(0)<<endl;
    Udm_.component(1) = pow(2,0.5)*pow((sigma_*g_.component(1)*(rhod_-rhoc_))/(pow(rhoc_,2)),0.25)*pow(1-alphad_,1.75);
    Info<<Udm_.component(1)<<endl;
    Udm_.component(2) = pow(2,0.5)*pow((sigma_*g_.component(2)*(rhod_-rhoc_))/(pow(rhoc_,2)),0.25)*pow(1-alphad_,1.75);
    Info<<Udm_.component(2)<<endl;
    Info<<((sigma_*g_*(rhoc_-rhod_))/(pow(rhoc_,2))).component(0)<<endl;*/
    
    Ddm_=(1-Mem_+VSMALL)*(rhoc_/rho())*D_/mag(alphad_);
    Ddm_.correctBoundaryConditions();
}


// ************************************************************************* //
