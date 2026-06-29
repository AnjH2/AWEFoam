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

#include "poreFlow.H"
#include "addToRunTimeSelectionTable.H"
#include "../../porousProperties/porousProperties.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace relativeVelocityModels
{
    defineTypeNameAndDebug(poreFlow, 0);
    addToRunTimeSelectionTable(relativeVelocityModel, poreFlow, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::relativeVelocityModels::poreFlow::poreFlow
(
    const dictionary& dict,
    const incompressibleTwoPhaseInteractingMixture& mixture,
    const word& modelName
)
:
    relativeVelocityModel(dict, mixture,modelName),
    continuousModel_(relativeVelocityModel::New(dict.subDict("contModel"),mixture_,modelName+"cont:")),
    dispersedModel_(relativeVelocityModel::New(dict.subDict("dispModel"),mixture_,modelName+"disp:")),
    electrodes({"Ne","Pe"}),
    Ddm_
    (
        IOobject
        (
            modelName+"Ddm",
            mixture_.alpha1().time().timeName(),
            mixture_.alpha1().mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mixture_.alpha1().mesh(),
        dimensionedTensor("Ddm_avg",dimVelocity*dimLength, Zero)
    ),
    alphaEff_
    (
        mixture_.alpha1().mesh().lookupObject<volScalarField>
        (
            "alphaEff"
        )
    ) 
{

}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::relativeVelocityModels::poreFlow::~poreFlow()
{}

// * * * * * * * * * * * * * * Private Functions  * * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
Foam::volScalarField Foam::relativeVelocityModels::poreFlow::chi()
{
    return pow(alphaEff_,2)*(3-2*alphaEff_); //Zero slope at both thresholds;
    //return pow(alphaEff_,3)*(10-15*alphaEff_+6*pow(alphaEff_,2)); // Zero first and second derivatives at both thresholds
}
void Foam::relativeVelocityModels::poreFlow::correct()
{    
   
    volScalarField chi_(chi());
    surfaceScalarField chif_(fvc::interpolate(chi_));
    dModel_->correct();
    continuousModel_->correct();
    dispersedModel_->correct();
    Info<<"after corrects()"<<endl;
    Udm_= (1-chi_)*continuousModel_->Udm()+chi_*dispersedModel_->Udm();
    Ddm_= (1-chi_)*continuousModel_->Ddm()+chi_*dispersedModel_->Ddm();
    BSCap_=(1-chi_)*continuousModel_->BSCap()+chi_*dispersedModel_->BSCap();
    Info<<"after Udm and Ddm"<<endl;
    F_=(1-chif_)*continuousModel_->F()+chif_*dispersedModel_->F();
    Info<<"after F"<<endl;
    //Ddm_=alphaMoveing_*baseModel_->Ddm();
}


// ************************************************************************* //
