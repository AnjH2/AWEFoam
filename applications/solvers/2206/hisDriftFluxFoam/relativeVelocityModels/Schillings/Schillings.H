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
    Foam::Schillings

Description
    Schillings relative velocity model

SourceFiles
    Schillings.C

\*---------------------------------------------------------------------------*/

#ifndef Schillings_H
#define Schillings_H

#include "relativeVelocityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace relativeVelocityModels
{

/*---------------------------------------------------------------------------*\
                         Class Schillings Declaration
\*---------------------------------------------------------------------------*/

class Schillings
:
    public relativeVelocityModel
{
    // Private data
	const incompressibleTwoPhaseInteractingMixture& mixture_;
	//const dimensionedVector g_;
	const meshObjects::gravity& g_;
	
	const dimensionedScalar rd_;
	const dimensionedScalar n_;
	const dimensionedScalar& rhoc_;
	const volScalarField& rhod_;
	dimensionedScalar D_;
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
    volVectorField USaff();
    volScalarField UHdiff();
    volScalarField USdiff();
    volVectorField USmig();
    
public:

    //- Runtime type information
    TypeName("Schillings");


    // Constructors

        //- Construct from components
        Schillings
        (
            const dictionary& dict,
            const incompressibleTwoPhaseInteractingMixture& mixture
        );


    //- Destructor
    ~Schillings();


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
