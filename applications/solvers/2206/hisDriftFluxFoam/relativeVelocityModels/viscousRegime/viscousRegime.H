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
    Foam::viscousRegime

Description
    viscousRegime relative velocity model

SourceFiles
    viscousRegime.C

\*---------------------------------------------------------------------------*/

#ifndef viscousRegime_H
#define viscousRegime_H

#include "relativeVelocityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace relativeVelocityModels
{

/*---------------------------------------------------------------------------*\
                         Class viscousRegime Declaration
\*---------------------------------------------------------------------------*/

class viscousRegime
:
    public relativeVelocityModel
{
    // Private data
	const incompressibleTwoPhaseInteractingMixture& mixture_;
	//const dimensionedVector g_;
	const meshObjects::gravity& g_;
	const dimensionedScalar sigma_;
	const dimensionedScalar rd_;
	const dimensionedScalar& rhoc_;
	const volScalarField& rhod_;
	dimensionedScalar D_;
	const volScalarField& V1_;
	const volScalarField& eps_;
	
    // Private functions??
    
    volScalarField rd();
    volScalarField psi();
    volScalarField f();
public:

    //- Runtime type information
    TypeName("viscousRegime");


    // Constructors

        //- Construct from components
        viscousRegime
        (
            const dictionary& dict,
            const incompressibleTwoPhaseInteractingMixture& mixture
        );


    //- Destructor
    ~viscousRegime();


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
