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

#include "combined.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace relativeVelocityModels
{
    defineTypeNameAndDebug(combined, 0);
    addToRunTimeSelectionTable(relativeVelocityModel, combined, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::relativeVelocityModels::combined::combined
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
    
    rd_(electrodes.size()),
	//(
	//	"R_DB",
	//	dimensionSet ( 0, 1, 0, 0, 0, 0,0),
	//	dict
	//),
    n_("n",dimless,dict_),
    minAlphad_("minAlphad",dimless,dict),
    rhoc_(mixture.rhoc()),
    rhod_(mixture.rhod()),
    
    CW_(electrodes.size()),//("CW", dimless, dict),
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
    V_
    (
        IOobject
        (
            "VCell",
            alphac_.time().timeName(),
            alphac_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alphac_.mesh(),
        dimensionedScalar(dimVolume, Zero)
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

    eg_("eg",(-1*g_)/mag(g_))
{
forAll(electrodes,i)
	{
	rd_.set
    	(
        	i,
        	new dimensionedScalar("rd_"+electrodes[i], dimLength,dict_)
        );
        CW_.set
        (
        	i,
        	new dimensionedScalar("CW_"+electrodes[i], dimless,dict_)
        );
        }
forAll (alphac_.mesh().C(), celli)
	{
	V_[celli]=alphac_.mesh().V()[celli];
	}        
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::relativeVelocityModels::combined::~combined()
{}

// * * * * * * * * * * * * * * Private Functions  * * * * * * * * * * * * * * //
volScalarField Foam::relativeVelocityModels::combined::kappa()
{
	
	return 60*pow(alphad_,1);
}

volScalarField Foam::relativeVelocityModels::combined::beta()
{
	
	return 1/3*pow(alphad_,2)*(1+0.5*exp(8.8*alphad_));
}

volScalarField Foam::relativeVelocityModels::combined::f()
{
	
	return pow(1-min(alphad_,minAlphad_),n_);
}

dimensionedVector Foam::relativeVelocityModels::combined::vStokes(int j)
{
	
	return mag(g_)*sqr(2*rd_[j])/(18*mixture_.nuc())*eg_;
}



volVectorField Foam::relativeVelocityModels::combined::UStokes(int j)
{
	
	return f()*mag(vStokes(j))*(1-(Ne_+Pe_)*hF_)*eg_;
}

volVectorField Foam::relativeVelocityModels::combined::ULub(int j)
{

/*forAll (alphac_.mesh().C(), celli)
	{
	ULub_[celli]=sign(yNormal_[celli])*
		(
				alphad_[celli]*CW_[j].value()*2/(2*rd_[j].value())*pow(mag(UStokes(j)[celli]),2)*pow(2*rd_[j].value()/(2*yNormal_[celli]),2)*alphac_.mesh().V()[celli]
		)
		/(
				3*Foam::constant::mathematical::pi*mixture_.nuc()()[celli]*2*rd_[j].value()
		)
			*vector(0.0,1.0,0.0);
	}
	
	
	return ULub_;*/
	/*return  -1*sign(yNormal_)*
			(
				alphad_*CW_[j]*2/(2*rd_[j])*pow(mag(UStokes(j)),2)*pow(2*rd_[j]/(2*yNormal_),2)
			)
			/(
				3*Foam::constant::mathematical::pi*mixture_.nuc()*2*rd_[j]
			)
			*dimensionedVector(dimless,vector(0.0,1.0,0.0));
				//Force/(3*pi*mu_L*D_db)*/
ULub_=sign(yNormal_)*
		(
				alphad_*CW_[j]*2/(2*rd_[j])*pow(mag(UStokes(j)),2)*pow(2*rd_[j]/(2*yNormal_),2)*V_
		)
		/(
				3*Foam::constant::mathematical::pi*mixture_.nuc()*2*rd_[j]
		)
			*dimensionedVector(dimless,vector(0.0,1.0,0.0));
return ULub_;
}




// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::relativeVelocityModels::combined::correct()
{    
    dModel_->correct();
    Udm_ = (1-Mem_+VSMALL)*(rhoc_/rho())*((UStokes(0)+ULub(0))*(Ne_+NeC_)+(UStokes(1)+ULub(1))*(Pe_+PeC_));

    /*
    Udm_.component(0) = pow(2,0.5)*pow((sigma_*g_.component(0)*(rhod_-rhoc_))/(pow(rhoc_,2)),0.25)*pow(1-alphad_,1.75);
    Info<<Udm_.component(0)<<endl;
    Udm_.component(1) = pow(2,0.5)*pow((sigma_*g_.component(1)*(rhod_-rhoc_))/(pow(rhoc_,2)),0.25)*pow(1-alphad_,1.75);
    Info<<Udm_.component(1)<<endl;
    Udm_.component(2) = pow(2,0.5)*pow((sigma_*g_.component(2)*(rhod_-rhoc_))/(pow(rhoc_,2)),0.25)*pow(1-alphad_,1.75);
    Info<<Udm_.component(2)<<endl;
    Info<<((sigma_*g_*(rhoc_-rhod_))/(pow(rhoc_,2))).component(0)<<endl;*/
    
    //Ddm_=(rhoc_/rho())*(rd_[0]*mag(vStokes(0))*D_[0]*(Ne_*dF_+NeC_)+rd_[1]*mag(vStokes(1))*D_[1]*(Pe_*dF_+PeC_))*f()/alphad_;
    //Ddm_.correctBoundaryConditions();
}


// ************************************************************************* //
