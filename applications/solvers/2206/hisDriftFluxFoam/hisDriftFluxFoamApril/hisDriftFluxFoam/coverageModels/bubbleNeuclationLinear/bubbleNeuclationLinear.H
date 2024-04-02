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
    Foam::bubbleNeuclationLinear

Description
    bubbleNeuclationLinear coverage model

SourceFiles
    bubbleNeuclationLinear.C

\*---------------------------------------------------------------------------*/

#ifndef bubbleNeuclationLinear_H
#define bubbleNeuclationLinear_H

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

class bubbleNeuclationLinear
:
    public coverageModel
{
    // Private data

        //- volume fraction of gas -> read
        const volScalarField& Psi_gas_;	//mass transfer from bubble neucelation 
        const volScalarField& Pe_;
        const volScalarField& Ne_;
        const volScalarField& J_;        
        volScalarField 	nb_;		//number of bubbles detaching
        volScalarField 	d_b_;
        volScalarField 	V_b_;
        volScalarField 	tau_d_;	//neucleation frequency
        volScalarField 	h_sb_;	//section height removed from sphere
        volScalarField 	V_sb_;	//volume of section removed frin sphere
        volScalarField 	A_sb_;	//surface area covered by a single bubble
        
        const PtrList <dimensionedScalar>& as_;		//specific surface area
        dimensionedScalar alpha_b_;	//contact angle
        
        dimensionedScalar d_ba_;
        dimensionedScalar d_bb_;
        
        dimensionedScalar tau_da_;
        dimensionedScalar tau_db_;
        label PI;
        
        
public:

    //- Runtime type information
    TypeName("bubbleNeuclationLinear");


    // Constructors

        //- Construct from components
        bubbleNeuclationLinear
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    ~bubbleNeuclationLinear();


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
