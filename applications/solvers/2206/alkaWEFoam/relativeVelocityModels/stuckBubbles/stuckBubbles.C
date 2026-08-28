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
#include "../../porousProperties/porousProperties.H"

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
    baseModel_(relativeVelocityModel::New(dict.subDict("baseModel"),mixture_,modelName+"base:")),
    covModel_(coverageModel::New(dict.subDict("covModel"),mixture_.alpha2().mesh(),mixture_,"base:")),
    electrodes({"Ne","Pe"}),
    dict_(dict),
    eps_
    (
        mixture_.alpha1().mesh().lookupObject<volScalarField>
        (
            "eps"
        )
    ),
    U_
    (
        mixture_.alpha1().mesh().lookupObject<volVectorField>
        (
            "U"
        )
    ),
    n_
	(
		"thetaInvPowN",
		dimless,
		 dict.lookupOrDefault<scalar>("thetaInvPowN", 1.0)
	),
	c_(electrodes.size()),
	rd_(electrodes.size()),
	Ddm_
    (
        IOobject
        (
            "Ddm_avg",
            U_.time().timeName(),
            U_.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        U_.mesh(),
        dimensionedTensor(dimVelocity*dimLength, Zero)
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
    ),

        as_(
	Ddm_.mesh().lookupObject<porousProperties>
        (
            "porousProperties"
        ).as()
        )
{
forAll(electrodes,i)
	{
	rd_.set
    	(
        	i,
        	new dimensionedScalar("rd_"+electrodes[i], dimLength,dict_)
        );
     c_.set
        (
            i,
		    new dimensionedScalar("thetaScale_"+electrodes[i],dimless,dict_)
	    );
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::relativeVelocityModels::stuckBubbles::~stuckBubbles()
{}

// * * * * * * * * * * * * * * Private Functions  * * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
// Local phase mass fractions used to combine stationary attached bubbles
// with the moving detached-bubble population.
Foam::volScalarField Foam::relativeVelocityModels::stuckBubbles::gamma_L()
{
    return alphac_*rhoc_/rho();
}
Foam::volScalarField Foam::relativeVelocityModels::stuckBubbles::gamma_G()
{
    return max(alphad_,SMALL)*rhod_/rho();
}
Foam::volScalarField Foam::relativeVelocityModels::stuckBubbles::gamma_GA()
{
    return alphaStatic_*rhod_/rho();
}
Foam::volScalarField Foam::relativeVelocityModels::stuckBubbles::gamma_GD()
{
    return alphaMoveing_*rhod_/rho();
}
void Foam::relativeVelocityModels::stuckBubbles::correct()
{    
    // Update both dispersion closures, the detached-bubble drift closure, and
    // the coverage model before partitioning the total gas saturation.
    dModel_->correct();
    baseModel_->correct();
    covModel_->correct();
    
    
    // Convert electrode coverage to attached-bubble saturation and cap it by
    // the total local gas saturation.
    alphaStatic_=min((c_[0]*rd_[0]*as_[0]*Ne_+c_[1]*rd_[1]*as_[1]*Pe_)*pow(covModel_->theta(),n_),alphad_);
    alphaMoveing_=alphad_-alphaStatic_;
    alphaStatic_.correctBoundaryConditions();
    // Reconstruct the overall gas relative velocity from stationary attached
    // and moving detached sub-phases. See review notes regarding the denominator.
    volScalarField scaleing=1/(gamma_G()*(1-gamma_GD()));
    Udm_ = scaleing*(-1*gamma_L()*gamma_GA()*U_+gamma_GD()*baseModel_->Udm());
    
    Ddm_=scaleing*(gamma_L()*gamma_GA()*dModel_->D()+gamma_GD()*baseModel_->Ddm());

}


// ************************************************************************* //
