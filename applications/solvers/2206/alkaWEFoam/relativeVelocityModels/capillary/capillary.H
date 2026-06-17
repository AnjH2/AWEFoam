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
    Foam::capillary

Description
    capillary relative velocity model,
    simplified model of the DarcyFlux.
    assuming that only the capillary pressure is drive the two phase flow.
        because the gas which is filling the pores are diconnected, 
        thus turning gravity and mixture pressure contributions to 0.
    both mobilities are non-zero.

SourceFiles
    capillary.C

\*---------------------------------------------------------------------------*/

#ifndef capillary_H
#define capillary_H

#include "relativeVelocityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace relativeVelocityModels
{

/*---------------------------------------------------------------------------*\
                         Class capillary Declaration
\*---------------------------------------------------------------------------*/

class capillary
:
    public relativeVelocityModel
{
	const incompressibleTwoPhaseInteractingMixture& mixture_;
	List<word> electrodes;
	const dictionary& dict_;
	
	volVectorField Urel_;
	dimensionedScalar CBeta_;
	
	
	const meshObjects::gravity& g_;
	const volScalarField& eps_;
	const volScalarField& K_;
	const volScalarField& kr_;
	const volScalarField& p_;
	const volScalarField& Solid_;
	
    // Private functions??
    volScalarField Mc();
    volScalarField Md();
    volScalarField Mm();
    volScalarField rhoS();
    volScalarField wd();
    volScalarField wc();
    
    



public:

    //- Runtime type information
    TypeName("capillary");


    // Constructors

        //- Construct from components
        capillary
        (
            const dictionary& dict,
            const incompressibleTwoPhaseInteractingMixture& mixture,
            const word& modelName = ""
        );


    //- Destructor
    ~capillary();


    // Member Functions
        //- Update the relative velocity
        void correctRelativeVelocity();
        //- Update the diffusion velocity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace relativeVelocityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
