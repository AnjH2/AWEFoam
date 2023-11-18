/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2016 OpenFOAM Foundation
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

#include "bubbleNeuclationLinear.H"
#include "addToRunTimeSelectionTable.H"
#include "../../porousProperties/porousProperties.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace coverageModels
{
    defineTypeNameAndDebug(bubbleNeuclationLinear, 0);
    addToRunTimeSelectionTable(coverageModel, bubbleNeuclationLinear, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coverageModels::bubbleNeuclationLinear::bubbleNeuclationLinear
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    coverageModel(dict,mesh),
    Psi_gas_
    (
        mesh.lookupObject<volScalarField>
        (
            "Psi_gas"
        )
    ),
        Pe_
    (
        mesh.lookupObject<volScalarField>
        (
            "Pe"
        )
    ),
        Ne_
    (
        mesh.lookupObject<volScalarField>
        (
            "Ne"
        )
    ),
            J_
    (
        mesh.lookupObject<volScalarField>
        (
            "J"
        )
    ),
   nb_(
        	IOobject
        	(
            	"nb",
            	mesh.time().timeName(),
            	mesh,
            	IOobject::NO_READ,
            	IOobject::AUTO_WRITE
        	),
        mesh,
        dimensionedScalar("nb__",dimensionSet(0,-3,-1,0,0,0,0),Zero)
   	),
   	
    d_b_
	(
        	IOobject
        	(
            	"d_b",
            	mesh.time().timeName(),
            	mesh,
            	IOobject::NO_READ,
            	IOobject::AUTO_WRITE
        	),
        mesh,
        dimensionedScalar("d_b__",dimensionSet(0,1,0,0,0,0,0),Zero)
   	),
    V_b_
	(
        	IOobject
        	(
            	"V_b",
            	mesh.time().timeName(),
            	mesh,
            	IOobject::NO_READ,
            	IOobject::AUTO_WRITE
        	),
        mesh,
        dimensionedScalar("V_b__",dimensionSet(0,3,0,0,0,0,0),Zero)
   	),

    tau_d_(
        	IOobject
        	(
            	"tau_d",
            	mesh.time().timeName(),
            	mesh,
            	IOobject::NO_READ,
            	IOobject::AUTO_WRITE
        	),
        mesh,
        dimensionedScalar("tau_d__",dimensionSet(0,0,-1,0,0,0,0),Zero)
   	),
    h_sb_
	(
        	IOobject
        	(
            	"h_sb",
            	mesh.time().timeName(),
            	mesh,
            	IOobject::NO_READ,
            	IOobject::AUTO_WRITE
        	),
        mesh,
        dimensionedScalar("h_sb__",dimensionSet(0,1,0,0,0,0,0),Zero)
   	),
    V_sb_
	(
        	IOobject
        	(
            	"V_sb",
            	mesh.time().timeName(),
            	mesh,
            	IOobject::NO_READ,
            	IOobject::AUTO_WRITE
        	),
        mesh,
        dimensionedScalar("V_sb__",dimensionSet(0,3,0,0,0,0,0),Zero)
   	),
    A_sb_	
	(
        	IOobject
        	(
            	"A_sb",
            	mesh.time().timeName(),
            	mesh,
            	IOobject::NO_READ,
            	IOobject::AUTO_WRITE
        	),
        mesh,
        dimensionedScalar("A_sb__",dimensionSet(0,2,0,0,0,0,0),Zero)
   	),
   as_(
	mesh.lookupObject<porousProperties>
        (
            "porousProperties"
        ).as()
        ),
  
     alpha_b_
	(
		"alpha_b",
		dimless,
		dict
	),
   d_ba_
	(
		"d_ba",
		dimensionSet(0,2,0,0,0,-1,0),
		dict
	),
	d_bb_
	(
		"d_bb",
		dimless,
		dict
	),
	tau_da_
	(
		"tau_da",
		dimensionSet(0,2,0,0,0,-1,0),
		dict
	),
	tau_db_
	(
		"tau_db",
		dimless,
		dict
	),
	PI(mesh.boundaryMesh().findPatchID("Outlet1"))
{

d_b_=(d_ba_*(mag(J_)/as_[0])+d_bb_)*dimensionedScalar("d__",dimensionSet(0,1,0,0,0,0,0),1);

V_b_=4/3.*Foam::constant::mathematical::pi*pow(d_b_/2,3);

h_sb_=d_b_/2-pow(pow(tan(alpha_b_)*d_b_,2)/(4*(1+pow(tan(alpha_b_),2))),0.5);

V_sb_=1./3*Foam::constant::mathematical::pi*pow(h_sb_,2)*(3*d_b_/2-h_sb_);

A_sb_=2*Foam::constant::mathematical::pi*d_b_/2*h_sb_;

tau_d_=(d_ba_*(mag(J_)/as_[0])+d_bb_)*dimensionedScalar("tau_d__",dimensionSet(0,0,-1,0,0,0,0),1);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::coverageModels::bubbleNeuclationLinear::~bubbleNeuclationLinear()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


void Foam::coverageModels::bubbleNeuclationLinear::correct()
{

  d_b_=((d_ba_*((mag(J_)+mag(J_.oldTime()))/2/as_[0])+d_bb_)*dimensionedScalar("d__",dimensionSet(0,1,0,0,0,0,0),1)+d_b_.oldTime())/2;

  V_b_=4/3.*Foam::constant::mathematical::pi*pow(d_b_/2,3);

  h_sb_=d_b_/2-pow(pow(tan(alpha_b_)*d_b_,2)/(4*(1+pow(tan(alpha_b_),2))),0.5);

  V_sb_=1./3*Foam::constant::mathematical::pi*pow(h_sb_,2)*(3*d_b_/2-h_sb_);

  A_sb_=2*Foam::constant::mathematical::pi*d_b_/2*h_sb_;

  nb_==-1.*Psi_gas_/(V_b_-V_sb_);

  tau_d_=(tau_da_*(mag(J_)/as_[0])+tau_db_)*dimensionedScalar("tau_d__",dimensionSet(0,0,-1,0,0,0,0),1);

 // Info<<"nb="<<nb_<<endl;
  //Info<<"numerator="<<A_sb_*nb_*tau_d_<<endl;
  theta_ ==max(1.-(A_sb_*nb_/(as_[0]*tau_d_)),Zero);

  //theta_ ==max(min(1-((A_sb_*nb_*tau_d_/as_[0])*Ne_+(A_sb_*nb_*tau_d_/as_[1])*Pe_),Pe_+Ne_),theta_*0);


  
}


// ************************************************************************* // 
