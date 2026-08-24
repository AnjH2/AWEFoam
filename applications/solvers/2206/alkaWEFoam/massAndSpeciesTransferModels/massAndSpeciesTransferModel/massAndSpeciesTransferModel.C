/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2017 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "massAndSpeciesTransferModel.H"
#include "surfaceFields.H"
#include "fvc.H"
#include "../../reactionProperties/reactionProperties.H"
#include "../../porousProperties/porousProperties.H"
#include "../../speciesProperties/speciesProperties.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(massAndSpeciesTransferModel, 0);
    defineRunTimeSelectionTable(massAndSpeciesTransferModel, dictionary);
}

// * * * * * * * * * * * * * Private Member Functions   * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::massAndSpeciesTransferModel::massAndSpeciesTransferModel
(
    const dictionary& dict,
    const fvMesh& mesh,
    const incompressibleTwoPhaseInteractingMixture& mixture
)
:
	mixture_(mixture),
	massAndSpeciesTransferModelDict_(dict),
	species2({"H2","O2","H2O","OH"}),
	species1({"H2","O2","H2O"}),
	speciesE({"H2","O2"}),
	
	Psi_BV_(species2.size()),
    
    psi_(
	mesh.lookupObject<reactionProperties>
        (
            "reactionProperties"
        ).psi()
        ),
    MW_(
	mesh.lookupObject<speciesProperties>
        (
            "speciesProperties"
        ).MW()
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
    PeC_
    (
        mesh.lookupObject<volScalarField>
        (
            "PeC"
        )
    ),
    NeC_
    (
        mesh.lookupObject<volScalarField>
        (
            "NeC"
        )
    ),
    J_
    (
        mesh.lookupObject<volScalarField>
        (
            "J"
        )
    ),
    T_
    (
        mesh.lookupObject<volScalarField>
        (
            "T"
        )
    ),
    alpha_
    (
        mesh.lookupObject<volScalarField>
        (
            "alpha.gas"
        )
    ),
    epsilon_
    (
        mesh.lookupObject<volScalarField>
        (
            "eps"
        )
    ),
        as_(
	mesh.lookupObject<porousProperties>
        (
            "porousProperties"
        ).as()
        ),
        
    mDot_Wall_(species2.size()),
    mDotAlpha_Wall_(species2.size()),
    waterVapour_(dict.get<bool>("waterVapour")),
    
    p_water_
    (
    	IOobject
    	(
    	    "p_water",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
    	),
    	mesh,
    	dimensionedScalar("p_water",dimPressure,1)
    ),
    p_water_pure_
    (
    	IOobject
    	(
    	    "p_water_pure",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
    	),
    	mesh,
    	dimensionedScalar("p_water_pure",dimPressure,1)
    ),
    C_sat_(species2.size())
 
{
forAll(species2,i)
	{
        mDot_Wall_.set
    	(
        	i,
        	new volScalarField
        	(
        		IOobject
        		(
        			"mDot_Wall_"+species2[i],
        			mesh.time().timeName(),
                		mesh,
                		IOobject::NO_READ,
                		IOobject::AUTO_WRITE
            		),
            		mesh,
            		dimensionSet(1,-3,-1,0,0,0,0)
        	)
        );
        mDotAlpha_Wall_.set
    	(
        	i,
        	new volScalarField
        	(
        		IOobject
        		(
        			"mDotAlpha_Wall_"+species2[i],
        			mesh.time().timeName(),
                		mesh,
                		IOobject::NO_READ,
                		IOobject::NO_WRITE
            		),
            		mesh,
            		dimensionSet(1,-3,-1,0,0,0,0)
        	)
        );

	Psi_BV_.set
    	(
        	i,
        	new volScalarField
        	(
        		IOobject
        		(
        			"Psi_BV_"+species2[i],
        			mesh.time().timeName(),
                		mesh,
                		IOobject::NO_READ,
                		IOobject::AUTO_WRITE
            		),
            		mesh,
            		dimensionSet(0,-3,-1,0,1,0,0)
        	)
        );
        C_sat_.set
    	(
        	i,
        	new volScalarField
        	(
        		IOobject
        		(
        			"C_sat_"+species2[i],
        			mesh.time().timeName(),
                		mesh,
                		IOobject::NO_READ,
                		IOobject::NO_WRITE
            		),
            		mesh,
            		dimensionedScalar(dimensionSet(0,-3,0,0,1,0,0),1)
        	)
        );
	}
if (waterVapour_){
	Info<<"		waterVapour active"<<endl;

	} else {
	Info<<"		waterVapour inactive"<<endl;
	}
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::massAndSpeciesTransferModel> Foam::massAndSpeciesTransferModel::New
(
    const dictionary& dict,
    const fvMesh& mesh,
    const incompressibleTwoPhaseInteractingMixture& mixture
)
{
    const word modelType(dict.get<word>(typeName));

    Info<< "Selecting massAndSpeciesTransfer model " << modelType << endl;

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "massAndSpeciesTransfer model",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << abort(FatalIOError);
    }

    return
        autoPtr<massAndSpeciesTransferModel>
        (
            ctorPtr
            (
                dict.optionalSubDict(modelType + "Coeffs"),
                mesh,
                mixture
            )
        );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::massAndSpeciesTransferModel::~massAndSpeciesTransferModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

// Convert the local volumetric current density to species production
// using reaction stoichiometry and Faraday's law (cf. paper Eq. 13).
void Foam::massAndSpeciesTransferModel::correct_Psi_BV(int i)
{

			Psi_BV_[i]=-((psi_[i]*Ne_+psi_[i+4]*Pe_)/2)*(J_)/Foam::constant::physicoChemical::F;

	
}
// Evaluate the local water vapour pressure in the KOH solution and in
// pure water. The vapour pressure is later used for gas-phase H2O transfer.
void Foam::massAndSpeciesTransferModel::correct_waterPartialPressure(const volScalarField& m_KOH)
{

		dimensionedScalar	Td_("T",dimensionSet(0,0,0,1,0,0,0),1);
		dimensionedScalar	Pd_("P",dimPressure,1e5);
		p_water_=pow(10,-0.01508*m_KOH-0.00167788*pow(m_KOH,2)+2.25887e-5*pow(m_KOH,3)+
		(1-0.0012062*m_KOH+5.6024e-4*pow(m_KOH,2)-7.8228e-6*pow(m_KOH,3))*
		(35.4462-3343.93/T_*Td_-10.9*log10(T_/Td_)+0.0041645*T_/Td_))*Pd_;
		
		p_water_pure_=pow(10,(35.4462-3343.93/T_*Td_-10.9*log10(T_/Td_)+0.0041645*T_/Td_))*Pd_;
		
		p_water_.correctBoundaryConditions();
        p_water_pure_.correctBoundaryConditions();
	
}







// Sum the interphase mass source over H2, O2 and, when enabled, H2O vapour.
Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::massAndSpeciesTransferModel::mDot() 
{

    Pair<tmp<volScalarField>> sumMDot(this->mDot(0,0)[0]*0,this->mDot(0,0)[0]*0);
    
    forAll(species1,i){
    	if (i<=1 or (waterVapour_ and i==2)){
    		
    	
    		
    	
    		Pair<tmp<volScalarField>> mDot = this->mDot(i,0);
    		sumMDot[0] = sumMDot[0] + mDot[0];
    		sumMDot[1] = sumMDot[1] + mDot[1];
    	}
    }
    return sumMDot;
}

// Total wall gas source limited by electrochemical production. The final
// factor accounts for water vapour required to saturate the generated gas.
const Foam::volScalarField
Foam::massAndSpeciesTransferModel::mDotAlpha_Wall()
{

    return
    (
    	-1*min(
    	    (mDotAlpha_Wall_[0]+mDotAlpha_Wall_[1])*(1-pow(alpha_,5))
    	    ,Psi_BV_[0]*MW_[0]+Psi_BV_[1]*MW_[1]
    	    )
    	*((MW_[2]*(Ne_/MW_[0]+Pe_/MW_[1]))*((mixture_.p_num()/(mixture_.p_num()-p_water_))-1)+1)
    );
}


/*bool Foam::massAndSpeciesTransferModel::read()
{
    if (regIOobject::read())
    {
        return true;
    }

    return false;
}*/

// ************************************************************************* //
