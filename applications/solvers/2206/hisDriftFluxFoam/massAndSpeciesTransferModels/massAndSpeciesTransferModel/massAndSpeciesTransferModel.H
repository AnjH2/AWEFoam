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

Class
    Foam::massAndSpeciesTransferModel

Description

SourceFiles
    massAndSpeciesTransferModel.C

\*---------------------------------------------------------------------------*/

#ifndef massAndSpeciesTransferModel_H
#define massAndSpeciesTransferModel_H

#include "fvCFD.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "physicoChemicalConstants.H"
#include "incompressibleTwoPhaseInteractingMixture.H"
#include "sherwoodModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class massAndSpeciesTransferModel Declaration
\*---------------------------------------------------------------------------*/

class massAndSpeciesTransferModel
{
    // Private Member Functions
	
        //- No copy construct
        massAndSpeciesTransferModel(const massAndSpeciesTransferModel&) = delete;

        //- No copy assignment
        void operator=(const massAndSpeciesTransferModel&) = delete;


protected:

    // Protected data
        //- Mixture properties
        const incompressibleTwoPhaseInteractingMixture& mixture_;
	const dictionary massAndSpeciesTransferModelDict_;
	List<word> species2;
	List<word> species1;
	List<word> speciesE;
	
	PtrList <volScalarField> Psi_BV_;
	//variables related to butlerVolmer species transfer 
	const PtrList <dimensionedScalar>& z_;
	const PtrList <dimensionedScalar>& psi_;
	const PtrList <dimensionedScalar>& MW_;
	const PtrList <volScalarField>& D2_;
	const volScalarField& Pe_;
        const volScalarField& Ne_; 
        const volScalarField& PeC_;
        const volScalarField& NeC_; 
        const volScalarField& J_;
        const volScalarField& T_;
        const volScalarField& alpha_;
        const volScalarField& epsilon_;
        const PtrList <dimensionedScalar>& as_;
        const dimensionedScalar& t_OH_;
        const volVectorField& U_;
        //const dimensionedScalar& p_num_;
        //Mass transfer model.
        PtrList <volScalarField> mDot_Wall_;
        PtrList <volScalarField> mDotAlpha_Wall_;
        PtrList <volScalarField> mDotDiff_;
        //water vapourModel
        const bool waterVapour_;
        
        volScalarField	p_water_;
        volScalarField	p_water_pure_;
       	PtrList <volScalarField> C_sat_;
       	const volScalarField& V2_;
        




public:

    //- Runtime type information
    TypeName("massAndSpeciesTransferModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        massAndSpeciesTransferModel,
        dictionary,
        (const dictionary& dict,
         const fvMesh& mesh,
         const incompressibleTwoPhaseInteractingMixture& mixture),
        (dict, mesh, mixture)
    );


    // Constructors

        //- Construct from components
        massAndSpeciesTransferModel
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const incompressibleTwoPhaseInteractingMixture& mixture
        );


    // Selector
    static autoPtr<massAndSpeciesTransferModel> New
    (
        const dictionary& dict,
        const fvMesh& mesh,
        const incompressibleTwoPhaseInteractingMixture& mixture
    );


    //- Destructor
    virtual ~massAndSpeciesTransferModel();


    // Member Functions
    
        
        
    	const incompressibleTwoPhaseInteractingMixture& mixture() const{
    		return mixture_;
    	}
	const PtrList <volScalarField>& Psi_BV() const
        {
            return Psi_BV_;
        }
        const PtrList <volScalarField>& mDot_Wall() const
        {
            return mDot_Wall_;
        }
        const volScalarField& p_water() const
        {
        	return p_water_;
        }
        const volScalarField& p_water_pure() const
        {
        	return p_water_pure_;
        }
        
        const bool& waterVapour() const
        {
        	return waterVapour_;
        }
        
        const PtrList <volScalarField>& C_sat() const
        {
            return C_sat_;
        }
        
        void correct_Psi_BV(int i);
        
        //
	void correct_waterPartialPressure(const volScalarField& m_KOH);
	
	const volScalarField p_water(const volScalarField m_KOH_);
	
	
        
        //- Update the massAndSpeciesTransfer
        virtual void correct_mDot_wall(const int i,const PtrList<volScalarField>& C2_s, const volScalarField& theta) = 0;
        
                //- Return the mass condensation and vaporisation rates as a
        //  coefficient to multiply (1 - alphal) for the condensation rate
        //  and a coefficient to multiply  alphal for the vaporisation rate
        //  of each species
        virtual Pair<tmp<volScalarField>> mDotAlphal(const int i) = 0;

        //- Return the mass condensation and vaporisation rates as coefficients
        virtual Pair<tmp<volScalarField>> mDot(const int i,const bool Write) = 0;

	Pair<tmp<volScalarField>> mDot();
	
	Pair<tmp<volScalarField>> mDotAlphal();
	
	const volScalarField mDotAlpha_Wall();
        
        //- Return the volumetric condensation and vaporisation rates as a
        //  coefficient to multiply (1 - alphal) for the condensation rate
        //  and a coefficient to multiply  alphal for the vaporisation rate
        //  of each species
        virtual Pair<tmp<volScalarField>> vDotAlphal();

        //- Return the volumetric condensation and vaporisation rates as
        //  coefficients of each species
        virtual Pair<tmp<volScalarField>> vDot();
        
        virtual Pair<tmp<volScalarField>> vDot_rhoC1();

        //- Correct the phaseChange model
        virtual void correct(const int i) = 0;

        //- Read the transportProperties dictionary and update
        //virtual bool read();
        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
