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
    Foam::surfaceSaturation

Description
    surfaceSaturation loss model, used in the butlerVolmer model

SourceFiles
    surfaceSaturation.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceSaturation_H
#define surfaceSaturation_H

#include "volFields.H"
#include "massAndSpeciesTransferModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace massAndSpeciesTransferModels
{

/*---------------------------------------------------------------------------*\
                         Class linear Declaration
\*---------------------------------------------------------------------------*/

class surfaceSaturation
:
    public massAndSpeciesTransferModel
{
    // Private data
	const dimensionedScalar	tau_c_;
	const PtrList<volScalarField>& C2_;
	const PtrList <volScalarField>& k_H_;
	const PtrList <dimensionedScalar>& MW_;
	
public:

    //- Runtime type information
    TypeName("surfaceSaturation");


    // Constructors

        //- Construct from components
        surfaceSaturation
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const incompressibleTwoPhaseInteractingMixture& mixture
        );


    //- Destructor
    ~surfaceSaturation();


    // Member Functions
        //- Update coverage
        virtual void correct_Psi_m( const int i, const PtrList<volScalarField>& C2_s);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace massAndSpeciesTransfers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
