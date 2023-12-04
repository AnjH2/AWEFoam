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

#include "surfaceSaturation.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "../../reactionProperties/reactionProperties.H"
#include "../../porousProperties/porousProperties.H"
#include "../../speciesProperties/speciesProperties.H"
#include "../../speciesTransport/speciesTransport.H"
#include "../../incompressibleTwoPhaseInteractingMixture/incompressibleTwoPhaseInteractingMixture.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace massAndSpeciesTransferModels
{
    defineTypeNameAndDebug(surfaceSaturation, 0);
    addToRunTimeSelectionTable(massAndSpeciesTransferModel, surfaceSaturation, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::massAndSpeciesTransferModels::surfaceSaturation::surfaceSaturation
(
    const dictionary& dict,
    const fvMesh& mesh,
    const incompressibleTwoPhaseInteractingMixture& mixture
)
:
    massAndSpeciesTransferModel(dict,mesh,mixture),
    tau_c_
	(
		"tau_c",
		dimensionSet ( 0, 0, 1, 0, 0, 0,0),
		dict
	),
	    tau_b_
	(
		"tau_b",
		dimensionSet ( 0, 0, 1, 0, 0, 0,0),
		dict
	),
	C2_(
	mesh.lookupObject<speciesProperties>
        (
            "speciesProperties"
        ).C2()
        ),
    k_H_(
	mesh.lookupObject<reactionProperties>
        (
            "reactionProperties"
        ).k_H()
        ),
    MW_(
	mesh.lookupObject<speciesProperties>
        (
            "speciesProperties"
        ).MW()
        )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::massAndSpeciesTransferModels::surfaceSaturation::~surfaceSaturation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::massAndSpeciesTransferModels::surfaceSaturation::correct_Psi_m(const int i, const PtrList<volScalarField>& C2_s)
{
	if (i<=1){
		Psi_m_[i]=(Pe_+Ne_) */*
		(
			
		max(min(
			(1)/tau_c_*MW_[i]*(C2_s[i]-k_H_[i]*(mixture_.p_num()-p_water_)),
			MW_[i]*(C2_[i]-k_H_[i]*(mixture_.p_num()-p_water_))*dimensionedScalar("__",dimensionSet ( 0, 0, -1, 0, 0, 0,0),1/C2_[i].mesh().time().deltaTValue())+Psi_BV_[i]*MW_[i]
			)
			,Psi_m_[i]*0) //(1-alpha1)*epsilon1**/
			(max((1)/tau_c_*MW_[i]*(C2_s[i]-k_H_[i]*(mixture_.p_num()))*(1-alpha_)*epsilon_,Psi_m_[i]*0)+min((1)/tau_b_*MW_[i]*(C2_[i]-k_H_[i]*(mixture_.p_num()))*alpha_*epsilon_,Psi_m_[i]*0));
		
	} else if (i==2 and waterVapour_) {
		Psi_m_[i]=((Psi_m_[0]/MW_[0]+Psi_m_[1]/MW_[1])*(mixture_.p_num()/(mixture_.p_num()-p_water_))-(Psi_m_[0]/MW_[0]+Psi_m_[1]/MW_[1]))*MW_[2];
		
	} else {
		Psi_m_[i]=Psi_m_[0]*0;
		Info<<"no Vapour Please"<<endl;
	}
	Psi_m_[i].correctBoundaryConditions();
}


// ************************************************************************* //
