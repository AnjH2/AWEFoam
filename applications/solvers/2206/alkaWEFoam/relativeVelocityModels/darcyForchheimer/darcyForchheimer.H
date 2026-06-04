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
    Foam::darcyForchheimer

Description
    darcyForchheimer relative velocity model,
    it defines properties for each electrode.
    lift only

SourceFiles
    darcyForchheimer.C

\*---------------------------------------------------------------------------*/

#ifndef darcyForchheimer_H
#define darcyForchheimer_H

#include "relativeVelocityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace relativeVelocityModels
{

/*---------------------------------------------------------------------------*\
                         Class darcyForchheimer Declaration
\*---------------------------------------------------------------------------*/

class darcyForchheimer
:
    public relativeVelocityModel
{
	const incompressibleTwoPhaseInteractingMixture& mixture_;
	List<word> electrodes;
	const dictionary& dict_;
	//const dimensionedVector g_;
	const meshObjects::gravity& g_;
	
	const volScalarField& rb_;
	PtrList <dimensionedScalar> CF_;
	dimensionedScalar n_;
	const dimensionedScalar minAlphad_;
	const dimensionedScalar& rhoc_;
	const volScalarField& rhod_;
	const volScalarField& eps_;
	const volScalarField& p_;
	const volScalarField& K_;
	
	
    // Private functions??
    
    
    volScalarField f();
    volScalarField A();
    volScalarField B();
    volVectorField force();
    volScalarField krd();
    volScalarField krc();
    



public:

    //- Runtime type information
    TypeName("darcyForchheimer");


    // Constructors

        //- Construct from components
        darcyForchheimer
        (
            const dictionary& dict,
            const incompressibleTwoPhaseInteractingMixture& mixture,
            const word& modelName = ""
        );


    //- Destructor
    ~darcyForchheimer();


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
