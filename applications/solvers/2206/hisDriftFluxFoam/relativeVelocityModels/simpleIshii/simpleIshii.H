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
    Foam::simpleIshii

Description
    simpleIshii relative velocity model

SourceFiles
    simpleIshii.C

\*---------------------------------------------------------------------------*/

#ifndef simpleIshii_H
#define simpleIshii_H

#include "relativeVelocityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace relativeVelocityModels
{

/*---------------------------------------------------------------------------*\
                         Class simpleIshii Declaration
\*---------------------------------------------------------------------------*/

class simpleIshii
:
    public relativeVelocityModel
{
	const incompressibleTwoPhaseInteractingMixture& mixture_;
	//const dimensionedVector g_;
	const meshObjects::gravity& g_;
	
	const dimensionedScalar rd_;
	const dimensionedScalar n_;
	const dimensionedScalar minAlphad_;
	const dimensionedScalar& rhoc_;
	const volScalarField& rhod_;
	dimensionedScalar D_;
	dimensionedScalar CW_;
	volVectorField ULub_;
	const volScalarField& V1_;
	const volScalarField& eps_;
	const volVectorField& U_;
	const dimensionedVector eg_;
	const dimensionedVector en_;
	const dimensionedVector em_;
	
	
    // Private functions??
    
    volScalarField kappa();
    volScalarField beta();
    volScalarField f();
    volVectorField vStokes();
    volScalarField gamma(const dimensionedVector i);
    
    volVectorField UStokes();
    
    volVectorField ULub();



public:

    //- Runtime type information
    TypeName("simpleIshii");


    // Constructors

        //- Construct from components
        simpleIshii
        (
            const dictionary& dict,
            const incompressibleTwoPhaseInteractingMixture& mixture
        );


    //- Destructor
    ~simpleIshii();


    // Member Functions

        //- Update the diffusion velocity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace relativeVelocityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
