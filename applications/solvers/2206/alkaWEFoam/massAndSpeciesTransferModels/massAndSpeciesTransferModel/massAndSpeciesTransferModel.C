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
    z_(
	mesh.lookupObject<speciesProperties>
        (
            "speciesProperties"
        ).z()
        ),
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
    D2_(
    	mesh.lookupObject<speciesProperties>
        	(
            		"speciesProperties"
        	).D2()
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
    t_OH_
    (
	mesh.lookupObject<reactionProperties>
        (
            "reactionProperties"
        ).t_OH()
        ),
    U_
    (
        mesh.lookupObject<volVectorField>
        (
            "U"
        )
    ),
    rb_
    (
        mesh.lookupObject<volScalarField>
        (
            "rb"
        )
    ),
    /*p_num_(
	mesh.lookupObject<incompressibleTwoPhaseInteractingMixture>
        (
            "transportProperties"
        ).p_num()
        ),    */
        
    mDot_Wall_(species2.size()),
    mDotAlpha_Wall_(species2.size()),
    mDotAlpha_FIELD_(species2.size()),
    waterVapour_(dict.get<bool>("waterVapour")),
    
    p_water_
    (
    	IOobject
    	(
    	    "p_water",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
    	),
    	mesh,
    	dimPressure
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
    	dimPressure
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
                		IOobject::READ_IF_PRESENT,
                		IOobject::NO_WRITE
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
        mDotAlpha_FIELD_.set
    	(
        	i,
        	new volScalarField
        	(
        		IOobject
        		(
        			"vDot_"+species2[i],
        			mesh.time().timeName(),
                		mesh,
                		IOobject::READ_IF_PRESENT,
                		IOobject::NO_WRITE
            		),
            		mesh,
            		//dimensionSet(1,-3,-1,0,0,0,0) //kg/m^3/s
            		dimensionSet(0,0,-1,0,0,0,0) // 1/s
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
                		IOobject::READ_IF_PRESENT,
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


void Foam::massAndSpeciesTransferModel::correct_Psi_BV(int i)
{

		//if(z_[i].value() == 0) //not charged species C_[i]
		//{
			Psi_BV_[i]=-((psi_[i]*Ne_+psi_[i+4]*Pe_)/2)*(J_)/Foam::constant::physicoChemical::F;
		//}
		//else //charged species C_[i]
		//{
		//	Psi_BV_[i]=-(t_OH_/z_[i].value()+(psi_[i]*Ne_+psi_[i+4]*Pe_)/2)*(J_)/Foam::constant::physicoChemical::F;
		//}
	
}
//calculates the partial pressure of water in a salt mixture
void Foam::massAndSpeciesTransferModel::correct_waterPartialPressure(const volScalarField& m_KOH)
{
	
	if (waterVapour_) {
	
		dimensionedScalar	Td_("T",dimensionSet(0,0,0,1,0,0,0),1);
		dimensionedScalar	Pd_("P",dimPressure,1e5);
		p_water_=pow(10,-0.01508*m_KOH-0.00167788*pow(m_KOH,2)+2.25887e-5*pow(m_KOH,3)+
		(1-0.0012062*m_KOH+5.6024e-4*pow(m_KOH,2)-7.8228e-6*pow(m_KOH,3))*
		(35.4462-3343.93/T_*Td_-10.9*log10(T_/Td_)+0.0041645*T_/Td_))*Pd_;
		
		p_water_pure_=pow(10,(35.4462-3343.93/T_*Td_-10.9*log10(T_/Td_)+0.0041645*T_/Td_))*Pd_;
		
		p_water_.correctBoundaryConditions();
		
		
	} else {
	
		p_water_=dimensionedScalar(dimPressure,0);
		p_water_pure_=dimensionedScalar(dimPressure,0);
	}
	
}

Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::massAndSpeciesTransferModel::vDotAlphal() 
{

    Pair<tmp<volScalarField>> sumVDotAlpha(this->mDotAlphal(0)[0]*dimensionedScalar(dimless/dimDensity,Zero),this->mDotAlphal(0)[0]*dimensionedScalar(dimless/dimDensity,Zero));

    forAll(species1,i){
	if (i<=1 or (waterVapour_ and i==2)){
		//const volScalarField rho1i(MW_[i]*(mixture_.p_num())/(Foam::constant::physicoChemical::R*T_));
	
    		volScalarField alphalCoeff
    		(
        		
        		-(
        			1.0/mixture_.rhod() + mixture_.alpha1()
       				*(1.0/mixture_.rhoc() - 1.0/mixture_.rhod())
       				//*(1/mixture_.rho())
       			)
       			
    		);
     		Pair<tmp<volScalarField>> mDotAlphal= this->mDotAlphal(i);
    		sumVDotAlpha[0] = sumVDotAlpha[0]+alphalCoeff*mDotAlphal[0];
    		sumVDotAlpha[1] = sumVDotAlpha[1]+alphalCoeff*mDotAlphal[1];
    	}
    }

    return sumVDotAlpha;
}


Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::massAndSpeciesTransferModel::vDot() 
{

    Pair<tmp<volScalarField>> sumVDot(this->mDot(0,0)[0]*dimensionedScalar(dimless/dimDensity,Zero),this->mDot(0,0)[0]*dimensionedScalar(dimless/dimDensity,Zero));;
    
        	volScalarField limitedAlpha1
    	(
        	min(max(mixture_.alpha1(), scalar(0)), scalar(1))
    	);
    	
    	volScalarField limitedAlpha2
    	(
        	min(max(mixture_.alpha2(), scalar(0)), scalar(1))
    	);
    
    forAll(species1,i){
    	if (i<=1 or (waterVapour_ and i==2)){
    		//const volScalarField rho1i(MW_[i]*(mixture_.p_num())/(Foam::constant::physicoChemical::R*T_));
    	
    		volScalarField pCoeff(1.0/mixture_.rhoc() - 1.0/mixture_.rhod());
    	
    		//volScalarField pCoeff(-1.0/mixture_.rho());
    	
    		Pair<tmp<volScalarField>> mDotAlphal = this->mDotAlphal(i);
    		sumVDot[0] = sumVDot[0] + limitedAlpha1*pCoeff*mDotAlphal[0];
    		sumVDot[1] = sumVDot[1] + limitedAlpha2*pCoeff*mDotAlphal[1];
    	}
    }
    return sumVDot;
}

Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::massAndSpeciesTransferModel::vDot_rhoC1() 
{

    Pair<tmp<volScalarField>> sumVDot(this->mDot(0,0)[0]*dimensionedScalar("zeroCoeff1",dimless/dimDensity,0.0),this->mDot(0,0)[0]*dimensionedScalar("zeroCoeff2",dimless/dimDensity,0.0));
    
    forAll(species1,i)
    {
    //Info<<i<<endl;
    	if (i<=1 or (waterVapour_ and i==2))
    	{
    	//Info<<"start0"<<endl;
    	//Info<<"waterVapour_"<<waterVapour_<<endl;
    		//const volScalarField rho1i(MW_[i]*(mixture_.p_num())/(Foam::constant::physicoChemical::R*T_));
    	    volScalarField alphalCoeff(1.0/mixture_.rhod() - mixture_.alpha1()*(1.0/mixture_.rhod() - 1.0/mixture_.rhoc()));
    		volScalarField pCoeff("pCoeff",- 1.0/mixture_.rhod());
    	
    		Pair<tmp<volScalarField>> mDotAlphal = this->mDotAlphal(i);
    		mDotAlpha_FIELD_[i]=mixture_.alpha1()*pCoeff*mDotAlphal[1]();
    		sumVDot[0] = sumVDot[0] - alphalCoeff*mDotAlphal[0];
    		sumVDot[1] = sumVDot[1] - alphalCoeff*mDotAlphal[1];
    		
    	}
    	    //Info<<"sumVDot1_Ne:"<<sum((- 1.0/mixture_.rhod()*this->mDotAlphal(i)[1]()*Ne_)())<<endl;
    		//Info<<"sumVDot1_Pe:"<<sum((- 1.0/mixture_.rhod()*this->mDotAlphal(i)[1]()*Pe_)())<<endl;
    		Info<<"Done1"<<endl;
    }
    //Info<<"Done2"<<endl;
    return sumVDot;
}

Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::massAndSpeciesTransferModel::mDot() 
{

    Pair<tmp<volScalarField>> sumMDot(this->mDot(0,0)[0]*0,this->mDot(0,0)[0]*0);
    
    forAll(species1,i){
    	if (i<=1 or (waterVapour_ and i==2)){
    		//const volScalarField rho1i(MW_[i]*(mixture_.p_num())/(Foam::constant::physicoChemical::R*T_));
    	
    		
    	
    		Pair<tmp<volScalarField>> mDot = this->mDot(i,0);
    		sumMDot[0] = sumMDot[0] + mDot[0];
    		sumMDot[1] = sumMDot[1] + mDot[1];
    	}
    }
    return sumMDot;
}
Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::massAndSpeciesTransferModel::mDotAlphal() 
{

    Pair<tmp<volScalarField>> sumMDotAlphal(this->mDotAlphal(0)[0]*0,this->mDotAlphal(0)[0]*0);
    
    forAll(species1,i){
    	if (i<=1 or (waterVapour_ and i==2)){
    		//const volScalarField rho1i(MW_[i]*(mixture_.p_num())/(Foam::constant::physicoChemical::R*T_));
    	
    		
    	
    		Pair<tmp<volScalarField>> mDotAlphal = this->mDotAlphal(i);
    		sumMDotAlphal[0] = sumMDotAlphal[0] + mDotAlphal[0];
    		sumMDotAlphal[1] = sumMDotAlphal[1] + mDotAlphal[1];
    	}
    }
    return sumMDotAlphal;
}

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
