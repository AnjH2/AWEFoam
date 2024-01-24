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

#include "mixedSaturation.H"
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
    defineTypeNameAndDebug(mixedSaturation, 0);
    addToRunTimeSelectionTable(massAndSpeciesTransferModel, mixedSaturation, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::massAndSpeciesTransferModels::mixedSaturation::mixedSaturation
(
    const dictionary& dict,
    const fvMesh& mesh,
    const incompressibleTwoPhaseInteractingMixture& mixture
)
:
    massAndSpeciesTransferModel(dict,mesh,mixture),
    K_AB_
	(
		"K_AB",
		dimensionSet ( 0, 1, -1, 0, 0, 0,0),
		dict
	),
    c_AB_
	(
		"c_AB",
		dimless,
		dict
	),
   K_DB_
	(
		"K_DB",
		dimensionSet ( 0, 1, -1, 0, 0, 0,0),
		dict
	),
   R_DB_
	(
		"R_DB",
		dimensionSet ( 0, 1, 0, 0, 0, 0,0),
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
        )
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::massAndSpeciesTransferModels::mixedSaturation::~mixedSaturation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::massAndSpeciesTransferModels::mixedSaturation::correct_mDot_wall(const int i, const PtrList<volScalarField>& C2_s, const volScalarField& theta)
{
	if (i<=1){
		C_sat_[i]=k_H_[i]*(mixture_.p_num());
		C_sat_[i].correctBoundaryConditions();
		mDot_Wall_[i]=(Pe_+Ne_)*min(max(c_AB_*K_AB_*as_[i]*theta*MW_[i]*(C2_s[i]-C_sat_[i]),mDot_Wall_[i]*0),Psi_BV_[i]*MW_[i]);
		mDot_Wall_[i].correctBoundaryConditions();
	} else if (i==2 and waterVapour_) {
		mDot_Wall_[i]=((mDot_Wall_[0]/MW_[0]+mDot_Wall_[1]/MW_[1])*(mixture_.p_num()/(mixture_.p_num()-p_water_))-(mDot_Wall_[0]/MW_[0]+mDot_Wall_[1]/MW_[1]))*MW_[2];
		mDot_Wall_[i].correctBoundaryConditions();		
	} else {
		mDot_Wall_[i]=mDot_Wall_[0]*0;
		Info<<"no Vapour Please"<<endl;
	}
}


Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::massAndSpeciesTransferModels::mixedSaturation::mDotAlphal(const int i) const
{
    if (i<= 1){
    	const dimensionedScalar C0(dimensionSet(0,-3,0,0,1,0,0), Zero);
    
    	//const volScalarField rho1i(MW_[i]*(mixture_.p_num())/(Foam::constant::physicoChemical::R*T_));


    	return Pair<tmp<volScalarField>>
    	(
        	(Pe_+Ne_)*K_DB_/R_DB_*epsilon_*alpha_*MW_[i]*max(C_sat_[i] - C2_[i], C0),
       		-(mDot_Wall_[i]+(Pe_+Ne_)*K_DB_/R_DB_*epsilon_*alpha_*MW_[i]*max(C2_[i] - C2_[i], C0))
    	);
    }
    else if (i==2 and waterVapour_) {
    	//Pair<tmp<volScalarField>> mDotAlphaH2 = this->mDotAlphal(0);
    	//Pair<tmp<volScalarField>> mDotAlphaO2 = this->mDotAlphal(1);
    	return Pair<tmp<volScalarField>>
    	(
        	((this->mDotAlphal(0)[0]/MW_[0]+this->mDotAlphal(1)[0]/MW_[1])*(mixture_.p_num()/(mixture_.p_num()-p_water_))-(this->mDotAlphal(0)[0]/MW_[0]+this->mDotAlphal(1)[0]/MW_[1]))*MW_[2],
       		-((this->mDotAlphal(0)[1]/MW_[0]+this->mDotAlphal(1)[1]/MW_[1])*(mixture_.p_num()/(mixture_.p_num()-p_water_))-(this->mDotAlphal(0)[1]/MW_[0]+this->mDotAlphal(1)[1]/MW_[1]))*MW_[2]
    	);
    } else {
    
    Pair<tmp<volScalarField>> mDotAlphal = this->mDotAlphal(i);
    	return Pair<tmp<volScalarField>>
    	(
        	mDotAlphal[0]*0,
       		mDotAlphal[0]*0
    	);
    }
    
}


Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::massAndSpeciesTransferModels::mixedSaturation::mDot(const int i) const
{

    if (i<=1) {
    	volScalarField limitedAlpha1
    	(
        	min(max(mixture_.alpha1(), scalar(0)), scalar(1))
    	);


    	const dimensionedScalar C0(dimensionSet(0,-3,0,0,1,0,0), Zero);
    
    	//const volScalarField rho1i(MW_[i]*(mixture_.p_num())/(Foam::constant::physicoChemical::R*T_));

    	volScalarField mDotE
    	(
        	"mDotE_"+species2[i], mDot_Wall_[i]+(Pe_+Ne_)*K_DB_/R_DB_*epsilon_*MW_[i]*limitedAlpha1*max(C2_[i] - C_sat_[i], C0)
    	);
    	volScalarField mDotC
    	(
        	"mDotC_"+species2[i], (Pe_+Ne_)*K_DB_/R_DB_*epsilon_*MW_[i]*limitedAlpha1*max(C_sat_[i] - C2_[i], C0)
    	);

    	if (limitedAlpha1.mesh().time().outputTime())
    	{
        	mDotC.write();
        	mDotE.write();
    	}

    	return Pair<Foam::tmp<Foam::volScalarField>>
    	(
        	tmp<volScalarField>(new volScalarField(mDotC)),
        	tmp<volScalarField>(new volScalarField(-mDotE))
    	);
    } else if (i==2 and waterVapour_) {
    	//Pair<tmp<volScalarField>> mDotH2 = this->mDot(0);
    	//Pair<tmp<volScalarField>> mDotO2 = this->mDot(1);
    	
    	return Pair<tmp<volScalarField>>
    	(
        	((this->mDot(0)[0]/MW_[0]+this->mDot(1)[0]/MW_[1])*(mixture_.p_num()/(mixture_.p_num()-p_water_))-(this->mDot(0)[0]/MW_[0]+this->mDot(1)[0]/MW_[1]))*MW_[2],
       		-((this->mDot(0)[1]/MW_[0]+this->mDot(1)[1]/MW_[1])*(mixture_.p_num()/(mixture_.p_num()-p_water_))-(this->mDot(0)[1]/MW_[0]+this->mDot(1)[1]/MW_[1]))*MW_[2]
    	);
    } else {
    
    Pair<tmp<volScalarField>> mDotAlphal = this->mDotAlphal(i);
    	return Pair<tmp<volScalarField>>
    	(
        	mDotAlphal[0]*0,
       		mDotAlphal[0]*0
    	);
    }
}





void Foam::massAndSpeciesTransferModels::mixedSaturation::correct()
{
}

/*
bool Foam::massAndSpeciesTransferModels::constant::read()
{
    if (massAndSpeciesTransferModel::read())
    {
        subDict(type() + "Coeffs").readEntry("coeffC", coeffC_);
        subDict(type() + "Coeffs").readEntry("coeffE", coeffE_);

        return true;
    }

    return false;
}*/

// ************************************************************************* //
