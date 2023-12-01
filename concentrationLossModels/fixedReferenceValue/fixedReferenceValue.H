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
    Foam::fixedReferenceValue

Description
    fixedReferenceValue loss model, used in the butlerVolmer model

SourceFiles
    fixedReferenceValue.C

\*---------------------------------------------------------------------------*/

#ifndef fixedReferenceValue_H
#define fixedReferenceValue_H

#include "volFields.H"
#include "concentrationLossModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace concentrationLossModels
{

/*---------------------------------------------------------------------------*\
                         Class linear Declaration
\*---------------------------------------------------------------------------*/

class fixedReferenceValue
:
    public concentrationLossModel
{
    // Private data
	PtrList <dimensionedScalar> C2_Ref_;
	PtrList <dimensionedScalar> n_Ne_;
	PtrList <dimensionedScalar> n_Pe_;
	const volScalarField& C2_s_H2_;
	const volScalarField& C2_s_O2_;
	const volScalarField& C2_s_H2O_;
	const volScalarField& C2_s_OH_;
public:

    //- Runtime type information
    TypeName("fixedReferenceValue");


    // Constructors

        //- Construct from components
        fixedReferenceValue
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    ~fixedReferenceValue();


    // Member Functions
        //- Update coverage
        virtual void correct(const PtrList <volScalarField>& C2_s_);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace concentrationLosss
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
