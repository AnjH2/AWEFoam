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
    Foam::mixedSaturation

Description
    mixedSaturation loss model, used in the butlerVolmer model

SourceFiles
    mixedSaturation.C

\*---------------------------------------------------------------------------*/

#ifndef mixedSaturation_H
#define mixedSaturation_H

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

class mixedSaturation
:
    public massAndSpeciesTransferModel
{
    // Private data
	const dimensionedScalar	K_AB_;
	const dimensionedScalar	c_AB_;
	const dimensionedScalar	K_DB_;
	const dimensionedScalar	R_DB_;
	const PtrList<volScalarField>& C2_;
	const PtrList <volScalarField>& k_H_;
	const PtrList <dimensionedScalar>& MW_;
	
public:

    //- Runtime type information
    TypeName("mixedSaturation");


    // Constructors

        //- Construct from components
        mixedSaturation
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const incompressibleTwoPhaseInteractingMixture& mixture
        );


    //- Destructor
    ~mixedSaturation();


    // Member Functions
        //- Update coverage
        virtual void correct_Psi_m( const int i, const PtrList<volScalarField>& C2_s, const volScalarField& theta);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace massAndSpeciesTransfers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
