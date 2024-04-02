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
    Foam::simpleIshiiTEST

Description
    simpleIshiiTEST relative velocity model

SourceFiles
    simpleIshiiTEST.C

\*---------------------------------------------------------------------------*/

#ifndef simpleIshiiTEST_H
#define simpleIshiiTEST_H

#include "relativeVelocityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace relativeVelocityModels
{

/*---------------------------------------------------------------------------*\
                         Class simpleIshiiTEST Declaration
\*---------------------------------------------------------------------------*/

class simpleIshiiTEST
:
    public relativeVelocityModel
{
    // Private data

        //- A coefficient
        dimensionedScalar n_;

        //- Drift velocity
        dimensionedVector V0_;
	const volScalarField& V1_;
	const volScalarField& eps_;
	dimensionedScalar D_;


public:

    //- Runtime type information
    TypeName("simpleIshiiTEST");


    // Constructors

        //- Construct from components
        simpleIshiiTEST
        (
            const dictionary& dict,
            const incompressibleTwoPhaseInteractingMixture& mixture
        );


    //- Destructor
    ~simpleIshiiTEST();


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
