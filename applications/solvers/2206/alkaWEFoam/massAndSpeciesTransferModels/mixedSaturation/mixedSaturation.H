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
    Foam::mixedSaturation

Description
    Gas-liquid mass-transfer model that treats wall-attached and detached
    bubbles separately. H2/O2 transfer is driven by supersaturation relative
    to Henry-law equilibrium and evaluated using Sherwood-based coefficients.
    Optional water-vapour transfer assumes saturated gas.

SourceFiles
    mixedSaturation.C

\*---------------------------------------------------------------------------*/

#ifndef mixedSaturation_H
#define mixedSaturation_H

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

class mixedSaturation
:
    public massAndSpeciesTransferModel
{
    // Private data
    	autoPtr<sherwoodModel> shModelAB_;
    	autoPtr<sherwoodModel> shModelDB_;
	const dimensionedScalar	c_AB_;
	const PtrList<volScalarField>& C2_;
	const PtrList <volScalarField>& k_H_;
	
public:

    //- Runtime type information
    TypeName("mixedSaturation");


    // Constructors

        //- Construct from components
        mixedSaturation
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const incompressibleTwoPhaseInteractingMixture& mixture
        );


    //- Destructor
    ~mixedSaturation();


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
