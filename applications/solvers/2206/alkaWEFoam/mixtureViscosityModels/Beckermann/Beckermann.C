/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2017 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "Beckermann.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace mixtureViscosityModels
{
    defineTypeNameAndDebug(Beckermann, 0);

    addToRunTimeSelectionTable
    (
        mixtureViscosityModel,
        Beckermann,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mixtureViscosityModels::Beckermann::Beckermann
(
    const word& name,
    const dictionary& viscosityPropertiesSub1,
    const dictionary& viscosityPropertiesSub2,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const word modelName
)
:
    mixtureViscosityModel(name, viscosityPropertiesSub1,viscosityPropertiesSub2, U, phi),
    BeckermannCoeffsSub1_(viscosityPropertiesSub1.optionalSubDict(modelName + "Coeffs")),
    BeckermannCoeffsSub2_(viscosityPropertiesSub2.optionalSubDict(modelName + "Coeffs")),
    mud_ref_(species1.size()),
    T_muref_(species1.size()),
    n_mu_(species1.size()),
    mud_(species1.size()),
        a_(
            	IOobject
            	(
               		"a_",
                	U.mesh().time().timeName(),
                	U.mesh(),
                	IOobject::NO_READ,
                	IOobject::NO_WRITE
            	),
            	U.mesh(),
            	dimensionedScalar(dimless,0)
        
        ),
    mud_m_(
            	IOobject
            	(
               		"mud",
                	U.mesh().time().timeName(),
                	U.mesh(),
                	IOobject::NO_READ,
                	IOobject::NO_WRITE
            	),
            	U.mesh(),
            	dimensionedScalar(dimDynamicViscosity,0)
        
        ),
    rhoc_("rho", dimDensity, BeckermannCoeffsSub2_),
    kr_("kr", dimless, BeckermannCoeffsSub1_),
    //BeckermannViscosityExponent_("exponent", dimless, BeckermannCoeffs_),
    //muMax_("muMax", dimDynamicViscosity, BeckermannCoeffs_),
    alpha_
    (
        U.mesh().lookupObject<volScalarField>
        (
            IOobject::groupName
            (
                viscosityPropertiesSub1.getOrDefault<word>("alpha", "alpha"),
                viscosityPropertiesSub1.dictName()
            )
        )
    )
{
forAll(species1,i)
	{
        mud_ref_.set
    	(
        	i,
        	new dimensionedScalar("mu_ref_"+species1[i], dimDynamicViscosity,BeckermannCoeffsSub1_)
        );
        T_muref_.set
    	(
        	i,
        	new dimensionedScalar("T_muref_"+species1[i], dimTemperature,BeckermannCoeffsSub1_)
        );
        n_mu_.set
    	(
        	i,
        	new dimensionedScalar("n_mu_"+species1[i], dimless,BeckermannCoeffsSub1_)
        );
        mud_.set
    		(
        		i,
        		new volScalarField
        		(
        			IOobject
        			(
        				"mud_"+species1[i],
        				U.mesh().time().timeName(),
                			U.mesh(),
                			IOobject::NO_READ,
                			IOobject::NO_WRITE
            			),
            			mud_ref_[i]*pow(T_/T_muref_[i],n_mu_[i])
        		)
        	);
        
	}


forAll(species1,i)
	{
	a_=a_*0;
	
	forAll(species1,j)
		{
		
		phi1_.set
    		(
        		i*3+j,
        		new volScalarField
        		(
        			IOobject
        			(
        				"phi1_"+species1[i]+"_"+species1[j],
        				U.mesh().time().timeName(),
                			U.mesh(),
                			IOobject::NO_READ,
                			IOobject::NO_WRITE
            			),
            			pow(1+pow(mud_[i]/mud_[j],0.5)*pow(MW_[j]/MW_[i],0.25),2)
            			/(4/pow(2,0.5)*pow(1+MW_[i]/MW_[j],0.5))
        		)
        	);
		}
	}
	mud_m_=(Mem_/2+(NeC_+Ne_))*mud_[0]+(Mem_/2+(PeC_+Pe_))*mud_[1];


}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::mixtureViscosityModels::Beckermann::mu(const volScalarField& muc, const volScalarField& rhod) const
{

    return pow(pow(1-min(alpha_,0.999),kr_)/(muc/rhoc_)+pow(min(alpha_,0.999),kr_)/(mud_m_/rhod),-1)*(min(alpha_,0.999)*rhod+(1-min(alpha_,0.999))*rhoc_);
}

void Foam::mixtureViscosityModels::Beckermann::mud_m_correct(const volScalarField& p_water,const dimensionedScalar& p_num)
{

mud_m_=mud_m_*0;

forAll(species1,i)
	{
	
	a_=a_*0;
	forAll(species1,j)
		{
		
			phi1_[i*3+j]=pow(1+pow(mud_[i]/mud_[j],0.5)*pow(MW_[j]/MW_[i],0.25),2)
            			*(1/pow(8,0.5)*pow(1+MW_[i]/MW_[j],-0.5));

        		if (j==0 and j!=i) {
        			a_=a_+(NeC_+Ne_)*phi1_[i*3+j]*(p_num-p_water)/p_num;
        		} else if(j==1 and j!=i) {
        			a_=a_+(PeC_+Pe_)*phi1_[i*3+j]*(p_num-p_water)/p_num;
        		} else if (j!=i) { 
        			a_=a_+Mem_*phi1_[i*3+j]+(Pe_+Ne_)*phi1_[i*3+j]*p_water/p_num;
        		}
		}
		if (i==0) {
			mud_m_=mud_m_+((NeC_+Ne_)*mud_[i])/(1+1/((p_num-p_water)/p_num)*a_);
		} else if(i==1) {
			mud_m_=mud_m_+((PeC_+Pe_)*mud_[i])/(1+1/((p_num-p_water)/p_num)*a_);
		} else { 
			mud_m_=mud_m_+mud_[i]/(1+1/((p_water)/p_num)*a_);
		}
	//mud_m_=mud_m_+(x_(i)*mud_[i])/a_;
	}


}


bool Foam::mixtureViscosityModels::Beckermann::read
(
    const dictionary& viscosityPropertiesSub1
)
{
    mixtureViscosityModel::read(viscosityPropertiesSub1);

    BeckermannCoeffsSub1_ = viscosityPropertiesSub1.optionalSubDict(typeName + "Coeffs");

    //BeckermannCoeffsSub1_.readEntry("mu", mud_);
  //  BeckermannCoeffs_.readEntry("n", BeckermannViscosityExponent_);
   // BeckermannCoeffs_.readEntry("muMax", muMax_);

    return true;
}


// ************************************************************************* //
