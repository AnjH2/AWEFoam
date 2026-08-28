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
    Foam::linearVogt

Description
    Combined coverage closure taking the larger of a saturation-based
    coverage estimate and the current-density-based Vogt correlation.

SourceFiles
    linearVogt.C

\*---------------------------------------------------------------------------*/

#ifndef linearVogt_H
#define linearVogt_H

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

class linearVogt
:
    public coverageModel
{
    // Private data

        //- exchange current density -> read
        const volScalarField& J_;
        const volScalarField& T_;
        const dimensionedScalar J_lim_;
        const dimensionedScalar J_scale_;
        const dimensionedScalar T_ref_;
        
        const PtrList <dimensionedScalar>& as_;
        const volScalarField& Pe_;
        const volScalarField& Ne_; 
        //- volume fraction of gas -> read
        const volScalarField& alphad_;
        const dimensionedScalar n_;
        
        
        
public:

    //- Runtime type information
    TypeName("linearVogt");


    // Constructors

        //- Construct from components
        linearVogt
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const incompressibleTwoPhaseInteractingMixture& mixture,
        const word& modelName = ""
        );


    //- Destructor
    ~linearVogt();


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
