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
    Composite attached/detached bubble model. Attached bubbles are stationary,
    while detached bubbles use a selectable base relative-velocity model. The
    total drift and dispersion are reconstructed from the two gas sub-phases.

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
	List<word> electrodes;
	const dictionary& dict_;
	const volScalarField& eps_;
	const volVectorField& U_;
	const dimensionedScalar n_; //powerlaw that converts coverage to saturation
	PtrList <dimensionedScalar> c_; //linearScale that converts coverage to saturation
	PtrList <dimensionedScalar> rd_;
	//values to check and print
	volTensorField Ddm_;
	volScalarField alphaStatic_;
	volScalarField alphaMoveing_;
	const PtrList <dimensionedScalar>& as_;
	
    // Private functions
        volScalarField gamma_L();
    volScalarField gamma_G();
    volScalarField gamma_GA();
    volScalarField gamma_GD();
    



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
