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

#include "DahlkildU.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace relativeVelocityModels
{
    defineTypeNameAndDebug(DahlkildU, 0);
    addToRunTimeSelectionTable(relativeVelocityModel, DahlkildU, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::relativeVelocityModels::DahlkildU::DahlkildU
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
    
    rd_
	(
		"R_DB",
		dimensionSet ( 0, 1, 0, 0, 0, 0,0),
		dict
	),
    n_("n",dimless,dict),
    rhoc_(mixture.rhoc()),
    rhod_(mixture.rhod()),
    
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
    eg_("eg",(-1*g_)/mag(g_))
{

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::relativeVelocityModels::DahlkildU::~DahlkildU()
{}


// * * * * * * * * * * * * * * Private Functions  * * * * * * * * * * * * * * //



volScalarField Foam::relativeVelocityModels::DahlkildU::f()
{
	
	return alphac_/(1/(alphac_+SMALL));
}


dimensionedVector Foam::relativeVelocityModels::DahlkildU::vStokes()
{
	
	return g_.value()*(sqr(2*rd_)/(3*mixture_.nuc()))*dimensionedScalar(dimAcceleration,1);
}



volVectorField Foam::relativeVelocityModels::DahlkildU::UStokes()
{
	
	return alphad_*f()*mag(vStokes())*eg_*(1-Mem_+VSMALL);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::relativeVelocityModels::DahlkildU::correct()
{
    dModel_->correct();
    Udm_ = (UStokes())*((alphac_)*rhoc_/(rho()*eps_));    
    Udm_.correctBoundaryConditions();
    
}


// ************************************************************************* //
