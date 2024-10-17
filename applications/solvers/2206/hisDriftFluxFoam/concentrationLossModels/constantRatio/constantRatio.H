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
    Foam::constantRatio

Description
    constantRatio loss model, used in the butlerVolmer model

SourceFiles
    constantRatio.C

\*---------------------------------------------------------------------------*/

#ifndef constantRatio_H
#define constantRatio_H

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

class constantRatio
:
    public concentrationLossModel
{
    // Private data

public:

    //- Runtime type information
    TypeName("constantRatio");


    // Constructors

        //- Construct from components
        constantRatio
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const massAndSpeciesTransferModel& mSTaPtr
        );


    //- Destructor
    ~constantRatio();


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
