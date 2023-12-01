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
    Foam::bubbleNeuclation

Description
    bubbleNeuclation coverage model

SourceFiles
    bubbleNeuclation.C

\*---------------------------------------------------------------------------*/

#ifndef bubbleNeuclation_H
#define bubbleNeuclation_H

#include "volFields.H"
#include "coverageModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace coverageModels
{

/*---------------------------------------------------------------------------*\
                         Class linear Declaration
\*---------------------------------------------------------------------------*/

class bubbleNeuclation
:
    public coverageModel
{
    // Private data

        //- volume fraction of gas -> read
        const volScalarField& Psi_gas_;	//mass transfer from bubble neucelation 
        const volScalarField& Pe_;
        const volScalarField& Ne_;
        volScalarField	nb_;		//number of bubbles detaching
        
        const PtrList <dimensionedScalar>& as_;		//specific surface area
        dimensionedScalar d_b_;		//bubble diameter
        dimensionedScalar V_b_;		//spehere volume
        dimensionedScalar alpha_b_;	//contact angle
        dimensionedScalar tau_d_;	//bubble neucleation time.
        dimensionedScalar h_sb_;	//section height removed from sphere
        dimensionedScalar V_sb_;	//volume of section removed frin sphere
        dimensionedScalar A_sb_;	//surface area covered by a single bubble
        label PI;
        
        
public:

    //- Runtime type information
    TypeName("bubbleNeuclation");


    // Constructors

        //- Construct from components
        bubbleNeuclation
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    ~bubbleNeuclation();


    // Member Functions

        //- Update coverage
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace coverageModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
