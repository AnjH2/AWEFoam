/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2015 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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

Class
    Foam::incompressibleTwoPhaseInteractingMixture

Description
    A two-phase incompressible transportModel for interacting phases
    requiring the direct evaluation of the mixture viscosity,
    e.g. activated sludge or slurry.

SourceFiles
    incompressibleTwoPhaseInteractingMixture.C

\*---------------------------------------------------------------------------*/

#ifndef incompressibleTwoPhaseInteractingMixture_H
#define incompressibleTwoPhaseInteractingMixture_H

#include "compressibleTransportModel.H"
#include "incompressible/viscosityModels/viscosityModel/viscosityModel.H"
#include "mixtureViscosityModel.H"
#include "twoPhaseMixture.H"
//#include "IOdictionary.H"
#include "physicoChemicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
         Class incompressibleTwoPhaseInteractingMixture Declaration
\*---------------------------------------------------------------------------*/

class incompressibleTwoPhaseInteractingMixture
:
    public compressibleTransportModel,
    public twoPhaseMixture
{
protected:
	dictionary dict_;
	List<word> species2;
	List<word> species1;
	const PtrList <volScalarField>& C1_;
	const volScalarField&	T_;
	const PtrList <dimensionedScalar>& MW_;


	autoPtr<mixtureViscosityModel> muModel_;
        autoPtr<viscosityModel> nucModel_;
        
        
        dimensionedScalar p_num_;
        volScalarField rhod_;
        dimensionedScalar rhoc_;
        
        //- Optional diameter of the dispersed phase particles
        dimensionedScalar dd_;

        //- Optional maximum dispersed phase-fraction (e.g. packing limit)
        scalar alphaMax_;

        const volVectorField& U_;
        const surfaceScalarField& phi_;

        volScalarField mu_;
        volScalarField nu_;
        
        const volScalarField& prgh_;
        
        const volScalarField& Pe_;
        const volScalarField& Ne_;
        const volScalarField& Mem_;


public:

    TypeName("incompressibleTwoPhaseInteractingMixture");


    // Constructors

        //- Construct from components
        incompressibleTwoPhaseInteractingMixture
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
    	    const dictionary& dict
        );


    //- Destructor
    virtual ~incompressibleTwoPhaseInteractingMixture() = default;


    // Member Functions


        //- Return const-access to the mixture viscosityModel
        const mixtureViscosityModel& muModel() const
        {
            return *muModel_;
        }

        //- Return const-access to the continuous-phase viscosityModel
        const viscosityModel& nucModel() const
        {
            return *nucModel_;
        }

        //- Return const-access to the dispersed-phase density
        const volScalarField& rhod() const
        {
            return rhod_;
        }

        //- Return const-access to continuous-phase density
        const dimensionedScalar& rhoc() const
        {
            return rhoc_;
        };

        //- Return the diameter of the dispersed-phase particles
        const dimensionedScalar& dd() const
        {
            return dd_;
        }

        //- Optional maximum phase-fraction (e.g. packing limit)
        //  Defaults to 1
        scalar alphaMax() const
        {
            return alphaMax_;
        }

        //- Return const-access to the mixture velocity
        const volVectorField& U() const
        {
            return U_;
        }
        const dimensionedScalar& p_num() const
        {
            return p_num_;
        }

        //- Return the dynamic mixture viscosity
        tmp<volScalarField> mu() const
        {
            return mu_;
        }
        
             const PtrList<volScalarField>& mud()
        {
            return muModel_->mud();
        }
	
	virtual tmp<volScalarField> muc() const
	{
		return rhoc_*nucModel_->nu();
	} 

        //- Return the dynamic mixture viscosity for patch
        virtual tmp<scalarField> mu(const label patchi) const
        {
            return mu_.boundaryField()[patchi];
        }

        //- Return the mixture density
        volScalarField rho() const
        {
            
            return alpha1_*rhod_+ alpha2_*rhoc_;// 
        }
        
        

        //- Return the mixture density for patch
        virtual tmp<scalarField> rho(const label patchi) const
        {
            return
                alpha1_.boundaryField()[patchi]*rhod_
              + alpha2_.boundaryField()[patchi]*rhoc_.value();//.fcValue();
        }

        //- Return the mixture viscosity
        virtual tmp<volScalarField> nu() const
        {
            return mu_/rho();
        }


        //- Return the mixture viscosity for patch
        virtual tmp<scalarField> nu(const label patchi) const
        {
            return mu_.boundaryField()[patchi]/rho(patchi);
        }

        //- Correct the laminar viscosity
        virtual void correct(const volScalarField& p_water)
        {
            muModel_->mud_m_correct(p_water,p_num_);
            mu_ = muModel_->mu(rhoc_*nucModel_->nu(), rhod_);
        }

        virtual void correct_rhod(const volScalarField& p_water)
        {
        	rhod_=Pe_*(p_num_)/(Foam::constant::physicoChemical::R*T_)*((p_num_-p_water)*(MW_[1])+p_water*MW_[2])/p_num_
			+ Ne_*(p_num_)/(Foam::constant::physicoChemical::R*T_)*((p_num_-p_water)*(MW_[0])+p_water*MW_[2])/p_num_
			+ Mem_*(p_num_)/(Foam::constant::physicoChemical::R*T_)*MW_[2];
        }
        virtual void correct_nu()
        {
        	nu_=mu_/rho();
        }
        /*virtual void correct_rhoc()
        {
        	rhoc_=C2_[0]*MW_[0]+C2_[1]*MW_[1]+C2_[2]*MW_[2]+C2_[3]*MW_[3]+C2_[3]*MW_[4];
        }*/

        
        
//	virtual void correct();
        //- Read base transportProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
