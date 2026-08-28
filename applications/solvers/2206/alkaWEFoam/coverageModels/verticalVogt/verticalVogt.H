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
    Foam::verticalVogt

Description
    Velocity-corrected Vogt coverage model. The base current-density
    correlation is reduced by an empirical factor accounting for removal
    of attached bubbles by electrolyte flow.

SourceFiles
    verticalVogt.C

\*---------------------------------------------------------------------------*/

#ifndef verticalVogt_H
#define verticalVogt_H

#include "volFields.H"
#include "unitConversion.H"
#include "coverageModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace coverageModels
{

/*---------------------------------------------------------------------------*\
                         Class linear Declaration
\*---------------------------------------------------------------------------*/

class verticalVogt
:
    public coverageModel
{
    // Private data

        //- exchange current density -> read
        List<word> electrodes;
        const volScalarField& J_;
        const volScalarField& T_;
        const dimensionedScalar J_lim_;
        const dimensionedScalar J_scale_;
        const dimensionedScalar T_ref_;
        PtrList <dimensionedScalar> CD_;    //drag coefficient for bubbles
        
        
        
        
        // references
        const PtrList <dimensionedScalar>& as_;
        const volScalarField& Pe_;
        const volScalarField& Ne_;
        const volVectorField& U_;
        
        // functions
        volScalarField correction();
        
        
public:

    //- Runtime type information
    TypeName("verticalVogt");


    // Constructors

        //- Construct from components
        verticalVogt
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const incompressibleTwoPhaseInteractingMixture& mixture,
            const word& modelName = ""
        );


    //- Destructor
    ~verticalVogt();


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
