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
    baseModel_(relativeVelocityModel::New(dict.subDict("baseModel"),mixture_,modelName+"base:")),
    covModel_(coverageModel::New(dict.subDict("covModel"),mixture_.alpha2().mesh(),mixture_,"base:")),
    electrodes({"Ne","Pe"}),
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
	thetaC_
	    (
        IOobject
        (
            "thetaScale",
            eps_.time().timeName(),
            eps_.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        eps_.mesh(),
        dimensionedScalar("thetaScale",dimless, Zero)
    ),
	rb_
    (
        mixture_.alpha1().mesh().lookupObject<volScalarField>
        (
            "rb"
        )
    ),
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
        dimensionedTensor("Ddm_avg",dimVelocity*dimLength, Zero)
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
        dimensionedScalar("alphaStatic",dimless, Zero)
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
        dimensionedScalar("alphaMoveing",dimless, Zero)
    ),
    alphaResidual_
	(
		"alphaResidual",
		dimless,
		dict.lookupOrDefault<scalar>("alphaResidual", 0)
	),
        as_(
	Ddm_.mesh().lookupObject<porousProperties>
        (
            "porousProperties"
        ).as()
        )
{
Info<< "thetaC_.name()       = " << thetaC_.name() << nl;
Info<< "thetaC_.writeOpt()   = " << thetaC_.writeOpt() << nl;
Info<< "thetaC_.registered?  = "
    << U_.mesh().foundObject<volScalarField>("thetaScale") << nl;
Info<< "thetaC_.instance()   = " << thetaC_.instance() << nl;
    if (!thetaC_.headerOk())
    {
    Info<< "Trying manual write of " << thetaC_.name() << nl;
    thetaC_ =
        dimensionedScalar("thetaScale_Ne",dimless,dict_)*Ne_
      + dimensionedScalar("thetaScale_Pe",dimless,dict_)*Pe_
      + dimensionedScalar("thetaScale_min", dimless, SMALL);
        
    }
    else
    {
    Info<< "Trying manual write of " << thetaC_.name() << nl;
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::relativeVelocityModels::stuckBubbles::~stuckBubbles()
{}

// * * * * * * * * * * * * * * Private Functions  * * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
Foam::volScalarField Foam::relativeVelocityModels::stuckBubbles::gamma_L()
{
    return alphac_*rhoc_/rho();
}
Foam::volScalarField Foam::relativeVelocityModels::stuckBubbles::gamma_G()
{
    return alphad_*rhod_/rho();
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
    dModel_->correct();
    baseModel_->correct();
    covModel_->correct();
    
    
    //volScalarField alphaStatithetaC_(pow(covModel_->theta(),1/n_));
    alphaStatic_=min(thetaC_*rb_*as_*pow(covModel_->theta(),n_),alphad_);
    //volScalarField alphaMoveing_(max(alphad_-alphaStatithetaC_,Zero));
    alphaMoveing_=alphad_-alphaStatic_;
    //volScalarField scaleing=1/(gamma_G()*(1-gamma_GA()));
    /*volScalarField scaleing=1/(alphad_*(rhothetaC_*alphathetaC_+rhod_*alphaMoveing_));
    Info<<average(scaleing)<<endl;
    Udm_ = scaleing*(-1*alphathetaC_*rhothetaC_*alphaStatithetaC_*U_+rho()*alphaMoveing_*baseModel_->Udm());
    Info<<average(mag(Udm_))<<endl;
    Ddm_=scaleing*(alphathetaC_*rhothetaC_*alphaStatithetaC_*dModel_->D()+rho()*alphaMoveing_*baseModel_->Ddm());
    */
    volScalarField scaleing=1/(gamma_G()*(1-gamma_GD()));
    Udm_ = scaleing*(-1*gamma_L()*gamma_GA()*U_+gamma_GD()*baseModel_->Udm());
    
    Ddm_=scaleing*(gamma_L()*gamma_GA()*dModel_->D()+gamma_GD()*baseModel_->Ddm());
    
    F_=baseModel_->F();
    //Ddm_=alphaMoveing_*baseModel_->Ddm();
}


// ************************************************************************* //
