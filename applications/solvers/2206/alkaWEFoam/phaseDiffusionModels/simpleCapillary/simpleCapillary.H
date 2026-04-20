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
    Foam::simpleCapillary

Description
    simpleCapillary relative velocity model

SourceFiles
    simpleCapillary.C

\*---------------------------------------------------------------------------*/

#ifndef simpleCapillary_H
#define simpleCapillary_H

#include "phaseDiffusionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseDiffusionModels
{

/*---------------------------------------------------------------------------*\
                         Class simpleCapillary Declaration
\*---------------------------------------------------------------------------*/

class simpleCapillary
:
    public phaseDiffusionModel
{
    // Private data
    
    const incompressibleTwoPhaseInteractingMixture& mixture_;
	const volScalarField& eps_;
	const volScalarField& K_;
	const volScalarField& Solid_;
	
	
    volScalarField A();
    volScalarField krd();
    volScalarField krc();


public:

    //- Runtime type information
    TypeName("simpleCapillary");


    // Constructors

        //- Construct from components
        simpleCapillary
        (
            const dictionary& dict,
            const incompressibleTwoPhaseInteractingMixture& mixture,
            const word& modelName = ""
        );


    //- Destructor
    ~simpleCapillary();


    // Member Functions

        //- Update the diffusion velocity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace phaseDiffusionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
