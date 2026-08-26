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
    Foam::relativeVelocityModel

Description

SourceFiles
    relativeVelocityModel.C

\*---------------------------------------------------------------------------*/

#ifndef relativeVelocityModel_H
#define relativeVelocityModel_H

#include "fvCFD.H"
#include "dictionary.H"
#include "incompressibleTwoPhaseInteractingMixture.H"
#include "runTimeSelectionTables.H"

#include "phaseDiffusionModel.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class relativeVelocityModel Declaration
\*---------------------------------------------------------------------------*/

class relativeVelocityModel
{
    // Private Member Functions

        //- Return the list of patchFieldTypes for Udm derived from U
        wordList UdmPatchFieldTypes() const;

        //- No copy construct
        relativeVelocityModel(const relativeVelocityModel&) = delete;

        //- No copy assignment
        void operator=(const relativeVelocityModel&) = delete;


protected:

    // Protected data

        //- Mixture properties
        const incompressibleTwoPhaseInteractingMixture& mixture_;
        const dictionary& dict_;
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

        //- Dispersed diffusion velocity
        mutable volVectorField Udm_;
        
        
        const scalar y0_;
        const scalar y1_;
	const scalar y2_;
	const scalar y3_;
	volScalarField yNormal_;
	
	const scalar hF_;
	
	const volScalarField& Pe_;
	
	const volScalarField& Ne_;
	const volScalarField& PeC_;
	
	const volScalarField& NeC_;
	
	const volScalarField& Mem_;
	
	const volScalarField& eps_;
	
	const word modelName_;
	
	autoPtr<phaseDiffusionModel> dModel_;

public:

    //- Runtime type information
    TypeName("relativeVelocityModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        relativeVelocityModel,
        dictionary,
        (const dictionary& dict,
         const incompressibleTwoPhaseInteractingMixture& mixture,
         const word& modelName),
        (dict, mixture, modelName)
    );


    // Constructors

        //- Construct from components
        relativeVelocityModel
        (
            const dictionary& dict,
            const incompressibleTwoPhaseInteractingMixture& mixture,
            const word& modelName = ""
        );


    // Selector
    static autoPtr<relativeVelocityModel> New
    (
        const dictionary& dict,
        const incompressibleTwoPhaseInteractingMixture& mixture,
        const word& modelName = ""
    );


    //- Destructor
    virtual ~relativeVelocityModel();


    // Member Functions

        //- Mixture properties
        const incompressibleTwoPhaseInteractingMixture& mixture() const
        {
            return mixture_;
        }
        virtual bool runDiffusion()
        {
        	return dModel_->runDiffusion();
        }
        //- Return the mixture mean density
        tmp<volScalarField> rho() const;

        //- Return the diffusion velocity of the dispersed phase
        const volVectorField& Udm() const
        {
            return Udm_;
        }
        //- Return the relative velocity of the continues phase
	tmp<volVectorField> Ucm() const;
	
	//- Return the relative velocity of the dispersed phase, volumesentric
	tmp<volVectorField> Udj() const;
	
	//- Return the relative diffusion of the dispersed phase, volumesentric
	virtual const volTensorField& Ddm() const
        {
            return dModel_->D();
        }

        //- Return the stress tensor due to the phase transport
        tmp<volTensorField> tauDm() const;

        //- Update the diffusion velocity
        virtual void correct() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
