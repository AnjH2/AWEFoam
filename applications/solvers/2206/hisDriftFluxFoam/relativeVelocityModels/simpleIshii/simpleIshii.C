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
    minAlphad_("minAlphad",dimless,dict),
    rhoc_(mixture.rhoc()),
    rhod_(mixture.rhod()),
   
    D_("D", dimVelocity*dimLength, dict),
    
    CW_("CW", dimless, dict),
        ULub_
    (
        IOobject
        (
            "ULub",
            alphac_.time().timeName(),
            alphac_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alphac_.mesh(),
        dimensionedVector(dimVelocity, Zero)
    ),
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
	
	return 60*pow(alphad_,1);
}

volScalarField Foam::relativeVelocityModels::simpleIshii::beta()
{
	
	return 1/3*pow(alphad_,2)*(1+0.5*exp(8.8*alphad_));
}

volScalarField Foam::relativeVelocityModels::simpleIshii::f()
{
	
	return pow(1-min(alphad_,minAlphad_),n_);
}

volVectorField Foam::relativeVelocityModels::simpleIshii::vStokes()
{
	
	return mag(g_)*sqr(2*rd_)/(18*mixture_.nuc())*eg_;
}

volScalarField Foam::relativeVelocityModels::simpleIshii::gamma(const dimensionedVector i)
{
	
	
	return fvc::grad(U_&eg_)&i;
}

volVectorField Foam::relativeVelocityModels::simpleIshii::UStokes()
{
	
	return f()*mag(vStokes())*eg_*(1-Mem_+VSMALL);
}

volVectorField Foam::relativeVelocityModels::simpleIshii::ULub()
{

	return  sign(yNormal_)*sqrt
		(
			(
			alphad_*rhoc_*CW_*2/(2*rd_)*pow(mag(UStokes()),2)*pow(2*rd_/(2*yNormal_),2)*dimensionedScalar(dimLength,1.0)
			)
			/rhoc_
		)*dimensionedVector(dimless,vector(0.0,1.0,0.0));
}




// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::relativeVelocityModels::simpleIshii::correct()
{    
	ULub_=ULub();
    Udm_ = (1-Mem_+VSMALL)*(rhoc_/rho())*(UStokes()+ULub_);
    

    /*
    Udm_.component(0) = pow(2,0.5)*pow((sigma_*g_.component(0)*(rhod_-rhoc_))/(pow(rhoc_,2)),0.25)*pow(1-alphad_,1.75);
    Info<<Udm_.component(0)<<endl;
    Udm_.component(1) = pow(2,0.5)*pow((sigma_*g_.component(1)*(rhod_-rhoc_))/(pow(rhoc_,2)),0.25)*pow(1-alphad_,1.75);
    Info<<Udm_.component(1)<<endl;
    Udm_.component(2) = pow(2,0.5)*pow((sigma_*g_.component(2)*(rhod_-rhoc_))/(pow(rhoc_,2)),0.25)*pow(1-alphad_,1.75);
    Info<<Udm_.component(2)<<endl;
    Info<<((sigma_*g_*(rhoc_-rhod_))/(pow(rhoc_,2))).component(0)<<endl;*/
    
    Ddm_=(1-Mem_+VSMALL)*(rhoc_/rho())*D_/alphad_;
    Ddm_.correctBoundaryConditions();
}


// ************************************************************************* //
