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

#include "bubbleNeuclation.H"
#include "addToRunTimeSelectionTable.H"
#include "../../porousProperties/porousProperties.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace coverageModels
{
    defineTypeNameAndDebug(bubbleNeuclation, 0);
    addToRunTimeSelectionTable(coverageModel, bubbleNeuclation, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coverageModels::bubbleNeuclation::bubbleNeuclation
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
   	
    as_(
	mesh.lookupObject<porousProperties>
        (
            "porousProperties"
        ).as()
        ),
    d_b_
	(
		"d_b",
		dimensionSet (0,1,0,0,0,0,0),
		dict
	),
    V_b_
	(
		"V_b",
		dimensionSet (0,3,0,0,0,0,0),
		4/3.*Foam::constant::mathematical::pi*pow(d_b_.value()/2,3)
	),
    alpha_b_
	(
		"alpha_b",
		dimless,
		dict
	),
    tau_d_
		(
		"tau_d",
		dimensionSet (0,0,1,0,0,0,0),
		dict
	),
    h_sb_
	(
		"h_sb",
		dimensionSet (0,1,0,0,0,0,0),
		d_b_.value()/2-pow(pow(tan(alpha_b_.value())*d_b_.value(),2)/(4*(1+pow(tan(alpha_b_.value()),2))),0.5)
	),
    V_sb_
	(
		"V_sb",
		dimensionSet (0,3,0,0,0,0,0),
		1./3*Foam::constant::mathematical::pi*pow(h_sb_.value(),2)*(3*d_b_.value()/2-h_sb_.value())
	),
    A_sb_	
	(
		"A_sb",
		dimensionSet(0,2,0,0,0,0,0),
		2*Foam::constant::mathematical::pi*d_b_.value()/2*h_sb_.value()
	),
	PI(mesh.boundaryMesh().findPatchID("Outlet1"))
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::coverageModels::bubbleNeuclation::~bubbleNeuclation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


void Foam::coverageModels::bubbleNeuclation::correct()
{
  //Info<<"Psi_gas="<<Psi_gas_<<endl;
  //Info<<A_sb_<<endl;
  //Info<<tau_d_<<endl;
  //Info<<as_[0]<<endl;
  nb_==-1.*Psi_gas_/(V_b_-V_sb_);
  
 // Info<<"nb="<<nb_<<endl;
  //Info<<"numerator="<<A_sb_*nb_*tau_d_<<endl;
  theta_ ==max(1.-(A_sb_*nb_*tau_d_/as_[0]),Zero);
  //Info<<"theta="<<theta_<<endl;
  //theta_ ==max(min(1-((A_sb_*nb_*tau_d_/as_[0])*Ne_+(A_sb_*nb_*tau_d_/as_[1])*Pe_),Pe_+Ne_),theta_*0);
}


// ************************************************************************* // 
