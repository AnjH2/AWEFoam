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

#include "stuckBubbles.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace relativeVelocityModels
{
    defineTypeNameAndDebug(stuckBubbles, 0);
    addToRunTimeSelectionTable(relativeVelocityModel, stuckBubbles, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::relativeVelocityModels::stuckBubbles::stuckBubbles
(
    const dictionary& dict,
    const incompressibleTwoPhaseInteractingMixture& mixture,
    const word& modelName
)
:
    relativeVelocityModel(dict, mixture,modelName),
    mixture_(mixture),
    baseModel_(relativeVelocityModel::New(dict.subDict("baseModel"),mixture_,"base:")),
    covModel_(coverageModel::New(dict.subDict("covModel"),mixture_.alpha2().mesh(),"base:")),
    dict_(dict),
    alpha1_(mixture_.alpha1()),
    U_
    (
        mixture_.alpha1().mesh().lookupObject<volVectorField>
        (
            "U"
        )
    ),
    n_
	(
		"n",
		dimless,
		dict
	),
	Udm_P
    (
        IOobject
        (
            "Udm_P",
            U_.time().timeName(),
            U_.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        U_.mesh(),
        dimensionedVector(dimVelocity, Zero)
    ),
    alphaStatic_
    (
        IOobject
        (
            "alphaStatic",
            U_.time().timeName(),
            U_.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        U_.mesh(),
        dimensionedScalar(dimless, Zero)
    ),
        alphaMoveing_
    (
        IOobject
        (
            "alphaMoveing",
            U_.time().timeName(),
            U_.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        U_.mesh(),
        dimensionedScalar(dimless, Zero)
    )
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::relativeVelocityModels::stuckBubbles::~stuckBubbles()
{}

// * * * * * * * * * * * * * * Private Functions  * * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::relativeVelocityModels::stuckBubbles::correct()
{    
    dModel_->correct();
    baseModel_->correct();
    covModel_->correct();
    
    //volScalarField alphaStatic_(pow(covModel_->theta(),1/n_));
    alphaStatic_=pow(covModel_->theta(),1/n_);
    //volScalarField alphaMoveing_(max(alpha1_-alphaStatic_,Zero));
    alphaMoveing_=max(alpha1_-alphaStatic_,Zero);
    Udm_ = (-U_*alphaStatic_+alphaMoveing_*baseModel_->Udm())/(alphaMoveing_+alphaStatic_);
    Udm_P=baseModel_->Udm();
    
    
    //Ddm_=alphaMoveing_*baseModel_->Ddm();
}


// ************************************************************************* //
