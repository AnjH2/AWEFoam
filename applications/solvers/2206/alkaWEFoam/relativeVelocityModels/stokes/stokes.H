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
    Foam::stokes

Description
    Buoyancy-driven Stokes drift closure with electrode-specific bubble radii
    and a gas-fraction hindrance factor. Phase dispersion is handled separately
    by phaseDiffusionModel.

SourceFiles
    stokes.C

\*---------------------------------------------------------------------------*/

#ifndef stokes_H
#define stokes_H

#include "relativeVelocityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace relativeVelocityModels
{

/*---------------------------------------------------------------------------*\
                         Class stokes Declaration
\*---------------------------------------------------------------------------*/

class stokes
:
    public relativeVelocityModel
{
	
	List<word> electrodes;
	
	const meshObjects::gravity& g_;
	
	PtrList <dimensionedScalar> rd_;
	dimensionedScalar n_;
	const dimensionedVector eg_;
	
	
    // Private functions??
    
    
    volScalarField f();
    dimensionedVector vStokes(int j);
    volVectorField UStokes(int j);
    



public:

    //- Runtime type information
    TypeName("stokes");


    // Constructors

        //- Construct from components
        stokes
        (
            const dictionary& dict,
            const incompressibleTwoPhaseInteractingMixture& mixture,
            const word& modelName = ""
        );


    //- Destructor
    ~stokes();


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
