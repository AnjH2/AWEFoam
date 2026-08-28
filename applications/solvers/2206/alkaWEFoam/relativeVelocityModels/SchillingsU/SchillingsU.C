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

#include "SchillingsU.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace relativeVelocityModels
{
    defineTypeNameAndDebug(SchillingsU, 0);
    addToRunTimeSelectionTable(relativeVelocityModel, SchillingsU, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::relativeVelocityModels::SchillingsU::SchillingsU
(
    const dictionary& dict,
    const incompressibleTwoPhaseInteractingMixture& mixture,
    const word& modelName
)
:
    relativeVelocityModel(dict, mixture,modelName),
    electrodes({"Ne","Pe"}),
    dict_(dict),
    g_(meshObjects::gravity::New(mixture.U().time())),
    
    rd_(electrodes.size()),
    n_("n",dimless,dict),
    eg_("eg",(-1*g_)/mag(g_))
{
forAll(electrodes,i)
	{
	rd_.set
    	(
        	i,
        	new dimensionedScalar("rd_"+electrodes[i], dimLength,dict_)
        );
        }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::relativeVelocityModels::SchillingsU::~SchillingsU()
{}


// * * * * * * * * * * * * * * Private Functions  * * * * * * * * * * * * * * //

volScalarField Foam::relativeVelocityModels::SchillingsU::kappa()
{
	
	return 0.6*pow(alphad_,2);
}

volScalarField Foam::relativeVelocityModels::SchillingsU::beta()
{
	
	return (1/3)*pow(alphad_,2)*(1.0+0.5*exp(8.8*alphad_));
}

volScalarField Foam::relativeVelocityModels::SchillingsU::f()
{
	
	return pow(1.0-alphad_,n_);
}
// Unit direction for the Saffman lift contribution, based on U x curl(U).
volVectorField Foam::relativeVelocityModels::SchillingsU::n()
{
    volVectorField omega = fvc::curl(U_);
    volVectorField lamb = U_^omega;
    return (lamb)/(mag(lamb)+dimensionedScalar(dimLength/dimTime/dimTime,SMALL));
}

volVectorField Foam::relativeVelocityModels::SchillingsU::vStokes()
{
	
	return -g_.value()*(sqr(2*(rd_[0]*Ne_+rd_[1]*Pe_))/(18.0*mixture_.nuc()))*dimensionedScalar(dimAcceleration,1);
}

volScalarField Foam::relativeVelocityModels::SchillingsU::gamma()
{
	    volTensorField gU_=fvc::grad(U_);
   volTensorField E_=0.5*(gU_+gU_.T());
	return sqrt(2*( E_ && E_ ));
}

volScalarField Foam::relativeVelocityModels::SchillingsU::tau()
{

	return mixture_.mu()*gamma()+dimensionedScalar(dimPressure,SMALL);
}

volVectorField Foam::relativeVelocityModels::SchillingsU::UStokes()
{
	
	return f()*mag(vStokes())*eg_*(1-Mem_+VSMALL);
}

volVectorField Foam::relativeVelocityModels::SchillingsU::USaff()
{
	return -f()*sign(gamma())*mag(vStokes())*6.46/(6.0*Foam::constant::mathematical::pi)*sqrt((sqr(rd_[0]*Ne_+rd_[1]*Pe_)*mag(gamma()))/mixture_.nuc())*n()*(1-Mem_+VSMALL);
}

// Shear-induced migration driven by gradients in the local shear stress.
// This is distinct from the explicit saturation-gradient dispersion tensor.
volVectorField Foam::relativeVelocityModels::SchillingsU::USmig()
{
	return -1/alphad_*pow((rd_[0]*Ne_+rd_[1]*Pe_),2)*mag(gamma())*kappa()*fvc::grad(tau())/tau();
}
// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::relativeVelocityModels::SchillingsU::correct()
{

    // Hydrodynamic phase dispersion is updated independently through dModel_.
    // Udm_ therefore contains only the remaining drift/migration contributions.
    dModel_->correct();
    Udm_ = (UStokes()+USaff()+USmig())*((alphac_)*rhoc_/rho());    
    Udm_.correctBoundaryConditions();
    
    

}


// ************************************************************************* //
