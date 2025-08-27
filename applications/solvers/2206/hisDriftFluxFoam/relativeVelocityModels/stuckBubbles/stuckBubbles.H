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

Class
    Foam::stuckBubbles

Description
    stuckBubbles relative velocity model,
    it defines properties for each electrode.
    diffusion and lift only

SourceFiles
    stuckBubbles.C

\*---------------------------------------------------------------------------*/

#ifndef stuckBubbles_H
#define stuckBubbles_H

#include "relativeVelocityModel.H"
#include "coverageModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace relativeVelocityModels
{

/*---------------------------------------------------------------------------*\
                         Class stuckBubbles Declaration
\*---------------------------------------------------------------------------*/

class stuckBubbles
:
    public relativeVelocityModel
{
	const incompressibleTwoPhaseInteractingMixture& mixture_;
	autoPtr<relativeVelocityModel> baseModel_;
	autoPtr<coverageModel> covModel_;
	const dictionary& dict_;
	const volScalarField& alpha1_;
	const volVectorField& U_;
	const dimensionedScalar n_; //powerlaw that converts coverage to saturation
	const dimensionedScalar c_; //linearScale that converts coverage to saturation
	//values to check and print
	volTensorField Ddm_;
	volScalarField alphaStatic_;
	volScalarField alphaMoveing_;
	
	
    // Private functions??
    
    



public:

    //- Runtime type information
    TypeName("stuckBubbles");


    // Constructors

        //- Construct from components
        stuckBubbles
        (
            const dictionary& dict,
            const incompressibleTwoPhaseInteractingMixture& mixture,
            const word& modelName = ""
        );


    //- Destructor
    ~stuckBubbles();


    // Member Functions

        //- Update the diffusion velocity
        virtual void correct();
        
        	virtual const volTensorField& Ddm() const
        {
            return Ddm_;
        }
        virtual bool runDiffusion()
        {
        	return dModel_->runDiffusion() or baseModel_->runDiffusion();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace relativeVelocityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
