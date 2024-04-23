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

#include "viscousRegime.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace relativeVelocityModels
{
    defineTypeNameAndDebug(viscousRegime, 0);
    addToRunTimeSelectionTable(relativeVelocityModel, viscousRegime, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::relativeVelocityModels::viscousRegime::viscousRegime
(
    const dictionary& dict,
    const incompressibleTwoPhaseInteractingMixture& mixture
)
:
    relativeVelocityModel(dict, mixture),
    mixture_(mixture),
    g_(meshObjects::gravity::New(mixture.U().time())),
    sigma_("SurfaceTension", dimForce/dimLength, dict),
    rd_
	(
		"R_DB",
		dimensionSet ( 0, 1, 0, 0, 0, 0,0),
		dict
	),
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
        )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::relativeVelocityModels::viscousRegime::~viscousRegime()
{}


// * * * * * * * * * * * * * * Private Functions  * * * * * * * * * * * * * * //

volScalarField Foam::relativeVelocityModels::viscousRegime::rd()
{

	return rd_*pow(rhoc_*mag(g_)*(rhoc_-rhod_)/pow(mixture_.muc(),2),1.0/3.0);
}

volScalarField Foam::relativeVelocityModels::viscousRegime::psi()
{
	return 0.55*pow(pow(1+0.08*pow(rd(),3.0),4.0/7.0)-1,0.75);
}

volScalarField Foam::relativeVelocityModels::viscousRegime::f()
{
	return pow(alphac_,0.5)*mixture_.muc()/mixture_.mu();
}
// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::relativeVelocityModels::viscousRegime::correct()
{
    Udm_ = (Pe_+Ne_)*(rhoc_/rho())*(Pe_+Ne_)*10.8*pow(mixture_.muc()*mag(g_)*(rhoc_-rhod_)/(pow(rhoc_,2.0)),1.0/3.0)*g_/mag(g_)
    			*(pow(alphac_,1.5)*f())/rd()
    			*(pow(psi(),4.0/3.0)*(1+psi()))/(1+psi()*pow(f(),6.0/7.0));
    /*
    Udm_.component(0) = pow(2,0.5)*pow((sigma_*g_.component(0)*(rhod_-rhoc_))/(pow(rhoc_,2)),0.25)*pow(1-alphad_,1.75);
    Info<<Udm_.component(0)<<endl;
    Udm_.component(1) = pow(2,0.5)*pow((sigma_*g_.component(1)*(rhod_-rhoc_))/(pow(rhoc_,2)),0.25)*pow(1-alphad_,1.75);
    Info<<Udm_.component(1)<<endl;
    Udm_.component(2) = pow(2,0.5)*pow((sigma_*g_.component(2)*(rhod_-rhoc_))/(pow(rhoc_,2)),0.25)*pow(1-alphad_,1.75);
    Info<<Udm_.component(2)<<endl;
    Info<<((sigma_*g_*(rhoc_-rhod_))/(pow(rhoc_,2))).component(0)<<endl;*/
    
    Ddm_=(1-Mem_+VSMALL)*(rhoc_/rho())*eps_*D_/alphad_;
    Ddm_.correctBoundaryConditions();
}


// ************************************************************************* //
