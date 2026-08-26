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
    Foam::phaseDiffusionModel

Description

SourceFiles
    phaseDiffusionModel.C

\*---------------------------------------------------------------------------*/

#ifndef phaseDiffusionModel_H
#define phaseDiffusionModel_H

#include "fvCFD.H"
#include "dictionary.H"
#include "incompressibleTwoPhaseInteractingMixture.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class phaseDiffusionModel Declaration
\*---------------------------------------------------------------------------*/

class phaseDiffusionModel
{
    // Private Member Functions

        //- Return the list of patchFieldTypes for Udm derived from U
        wordList UdmPatchFieldTypes() const;

        //- No copy construct
        phaseDiffusionModel(const phaseDiffusionModel&) = delete;

        //- No copy assignment
        void operator=(const phaseDiffusionModel&) = delete;


protected:

    // Protected data

        //- Mixture properties
        const incompressibleTwoPhaseInteractingMixture& mixture_;

        //- Name of the continuous phase
        const word continuousPhaseName_;

        //- Continuous phase fraction
        const volScalarField& alphac_;

        //- Dispersed phase fraction
        const volScalarField& alphad_;

        //- Continuous density
        const dimensionedScalar& rhoc_;

        //- Dispersed density
        const volScalarField& rhod_;
        
        //- Dispersed phase diffusion        
        volTensorField Ddm_;
	
	const scalar dF_;
	
	const volScalarField& Pe_;
	
	const volScalarField& Ne_;
	const volScalarField& PeC_;
	
	const volScalarField& NeC_;
	
	const volScalarField& Mem_;
	
	const volScalarField& eps_;
	
	const word modelName_;
	
	volScalarField VToM();
	volScalarField RToM();
	
	

public:

    //- Runtime type information
    TypeName("phaseDiffusionModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        phaseDiffusionModel,
        dictionary,
        (const dictionary& dict,
         const incompressibleTwoPhaseInteractingMixture& mixture,
         const word& modelName),
        (dict, mixture, modelName)
    );


    // Constructors

        //- Construct from components
        phaseDiffusionModel
        (
            const dictionary& dict,
            const incompressibleTwoPhaseInteractingMixture& mixture,
            const word& modelName = ""
        );


    // Selector
    static autoPtr<phaseDiffusionModel> New
    (
        const dictionary& dict,
        const incompressibleTwoPhaseInteractingMixture& mixture,
        const word& modelName = ""
    );


    //- Destructor
    virtual ~phaseDiffusionModel();


    // Member Functions

        //- Mixture properties
        const incompressibleTwoPhaseInteractingMixture& mixture() const
        {
            return mixture_;
        }
	// returns if diffusion equations should be run, default is true,
	virtual bool runDiffusion()
	{
	    return true;
	}

	//- Return the relative diffusion of the dispersed phase, volumesentric
	const volTensorField& D() const
        {
            return Ddm_;
        }
	//- Return the relative diffusion of the continues phase
	tmp<volTensorField> Dcm() const;

        //- Update the diffusion velocity
        virtual void correct() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
