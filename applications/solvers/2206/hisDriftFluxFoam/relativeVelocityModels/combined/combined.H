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
    Foam::combined

Description
    combined relative velocity model,
    it defines properties for each electrode.

SourceFiles
    combined.C

\*---------------------------------------------------------------------------*/

#ifndef combined_H
#define combined_H

#include "relativeVelocityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace relativeVelocityModels
{

/*---------------------------------------------------------------------------*\
                         Class combined Declaration
\*---------------------------------------------------------------------------*/

class combined
:
    public relativeVelocityModel
{
	const incompressibleTwoPhaseInteractingMixture& mixture_;
	List<word> electrodes;
	const dictionary& dict_;
	//const dimensionedVector g_;
	const meshObjects::gravity& g_;
	
	PtrList <dimensionedScalar> rd_;
	dimensionedScalar n_;
	const dimensionedScalar minAlphad_;
	const dimensionedScalar& rhoc_;
	const volScalarField& rhod_;
	PtrList <dimensionedScalar> D_;
	PtrList <dimensionedScalar> CW_;
	volVectorField ULub_;
	const volScalarField& V1_;
	const volScalarField& eps_;
	const volVectorField& U_;
	const dimensionedVector eg_;
	
	
    // Private functions??
    
    volScalarField kappa();
    volScalarField beta();
    volScalarField f();
    volVectorField vStokes(int j);
    volVectorField UStokes(int j);
    
    volVectorField ULub(int j);



public:

    //- Runtime type information
    TypeName("combined");


    // Constructors

        //- Construct from components
        combined
        (
            const dictionary& dict,
            const incompressibleTwoPhaseInteractingMixture& mixture
        );


    //- Destructor
    ~combined();


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
