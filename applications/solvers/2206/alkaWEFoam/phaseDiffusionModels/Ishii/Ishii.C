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

\*---------------------------------------------------------------------------*/

#include "Ishii.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseDiffusionModels
{
    defineTypeNameAndDebug(Ishii, 0);
    addToRunTimeSelectionTable(phaseDiffusionModel, Ishii, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseDiffusionModels::Ishii::Ishii
(
    const dictionary& dict,
    const incompressibleTwoPhaseInteractingMixture& mixture,
    const word& modelName
)
:
    phaseDiffusionModel(dict, mixture,modelName),
    D_("DiffusionConstant", dimLength*dimLength/dimTime, dict),
    DT_(dimLength*dimLength/dimTime, Zero)
{
    Ddm_.replace(tensor::XX,D_);
    Ddm_.replace(tensor::YY,D_);
    Ddm_.replace(tensor::ZZ,D_);
    DT_.replace(tensor::XX,D_);
    DT_.replace(tensor::YY,D_);
    DT_.replace(tensor::ZZ,D_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseDiffusionModels::Ishii::~Ishii()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::phaseDiffusionModels::Ishii::correct()
{
    Ddm_=dF_*(1-Mem_)*VToM()*DT_/alphad_;
}


// ************************************************************************* //
