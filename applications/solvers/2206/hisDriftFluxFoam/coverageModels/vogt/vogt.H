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
    Foam::vogt

Description
    vogt coverage model

SourceFiles
    vogt.C

\*---------------------------------------------------------------------------*/

#ifndef vogt_H
#define vogt_H

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

class vogt
:
    public coverageModel
{
    // Private data

        //- exchange current density -> read
        const volScalarField& J_;
        const dimensionedScalar J_lim_;
        const dimensionedScalar J_scale_;
        
        const PtrList <dimensionedScalar>& as_;
        const volScalarField& Pe_;
        const volScalarField& Ne_;
        const volScalarField& p_water_;
        
        
        
public:

    //- Runtime type information
    TypeName("vogt");


    // Constructors

        //- Construct from components
        vogt
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    ~vogt();


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
