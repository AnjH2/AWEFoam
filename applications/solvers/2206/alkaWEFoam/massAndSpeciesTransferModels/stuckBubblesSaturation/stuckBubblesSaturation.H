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
    Foam::stuckBubblesSaturation

Description
    stuckBubblesSaturation mass transfer that have to be coupled with the
    relativ velocity model stuckBubbles.
    it considers both mass transfer to growing and moving bubbles seperatly.

SourceFiles
    stuckBubblesSaturation.C

\*---------------------------------------------------------------------------*/

#ifndef stuckBubblesSaturation_H
#define stuckBubblesSaturation_H

#include "volFields.H"
#include "massAndSpeciesTransferModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace massAndSpeciesTransferModels
{

/*---------------------------------------------------------------------------*\
                         Class linear Declaration
\*---------------------------------------------------------------------------*/

class stuckBubblesSaturation
:
    public massAndSpeciesTransferModel
{
    // Private data
    	autoPtr<sherwoodModel> shModelAB_;
    	autoPtr<sherwoodModel> shModelDB_;
	const dimensionedScalar	c_AB_;
	const PtrList<volScalarField>& C2_;
	const PtrList <volScalarField>& k_H_;
	const volScalarField& alphaStatic_;
	const volScalarField& alphaMoveing_;
	
public:

    //- Runtime type information
    TypeName("stuckBubblesSaturation");


    // Constructors

        //- Construct from components
        stuckBubblesSaturation
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const incompressibleTwoPhaseInteractingMixture& mixture
        );


    //- Destructor
    ~stuckBubblesSaturation();


    // Member Functions
        

        
        //- Update coverage        
        virtual void correct_mDot_wall( const int i, const PtrList<volScalarField>& C2_s, const volScalarField& theta);
        
        //- Return the mass condensation and vaporisation rates as a
        //  coefficient to multiply (1 - alphal) for the condensation rate
        //  and a coefficient to multiply  alphal for the vaporisation rate
        virtual Pair<tmp<volScalarField>> mDotAlphal(const int i);

        //- Return the mass condensation and vaporisation rates as coefficients
        virtual Pair<tmp<volScalarField>> mDot(const int i,const bool Write);

	
	
        //- Correct the constant phaseChange model
        virtual void correct(const int i);

        //- Read the transportProperties dictionary and update
        //virtual bool read();
        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace massAndSpeciesTransfers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
