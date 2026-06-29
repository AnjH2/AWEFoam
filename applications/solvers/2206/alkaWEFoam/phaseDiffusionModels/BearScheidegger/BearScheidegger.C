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

#include "BearScheidegger.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseDiffusionModels
{
    defineTypeNameAndDebug(BearScheidegger, 0);
    addToRunTimeSelectionTable(phaseDiffusionModel, BearScheidegger, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseDiffusionModels::BearScheidegger::BearScheidegger
(
    const dictionary& dict,
    const incompressibleTwoPhaseInteractingMixture& mixture,
    const word& modelName
)
:
    phaseDiffusionModel(dict, mixture,modelName),
    electrodes({"Ne","Pe"}),
    dict_(dict),
    g_(meshObjects::gravity::New(mixture.U().time())),
    U_(
	        alphad_.mesh().lookupObject<volVectorField>
        	(
            		"U"
        	)
	),
    //dimensionedScalar
    dispT_ 
    (
            IOobject
            (
                "dispT",
                alphad_.mesh().time().timeName(),
                alphad_.mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE         // run-time field, not written
            ),
            alphad_.mesh(),
            dimensionedScalar("dispT", dimLength,0 ) // start at 0
    ),
    dispL_ 
    (
            IOobject
            (
                "dispL",
                alphad_.mesh().time().timeName(),
                alphad_.mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE         // run-time field, not written
            ),
            alphad_.mesh(),
            dimensionedScalar("dispL", dimLength,0 ) // start at 0
    )
    
    
{
    if (!dispT_.headerOk())
    {
    dispT_ =
        dimensionedScalar("dispT_Ne",dimLength,dict_)*(Ne_+NeC_)
      + dimensionedScalar("dispT_Pe",dimLength,dict_)*(Pe_+PeC_);
        
    }
    if (!dispL_.headerOk())
    {
    dispL_ =
        dimensionedScalar("dispL_Ne",dimLength,dict_)*(Ne_+NeC_)
      + dimensionedScalar("dispL_Pe",dimLength,dict_)*(Pe_+PeC_);
        
    }
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseDiffusionModels::BearScheidegger::~BearScheidegger()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


volTensorField Foam::phaseDiffusionModels::BearScheidegger::UBSdiff()
{
const dimensionedScalar Umin ( "Umin", dimVelocity, SMALL );
volScalarField magU(mag(U_));
    return //Transverse isotropic contribution 
            dispT_*magU*tensor::I 
            // Additional longitudinal contribution 
            + (dispL_ - dispT_) *sqr(U_) /max(magU, Umin);
}



// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::phaseDiffusionModels::BearScheidegger::correct()
{

    Ddm_=(1-Mem_)*(VToM()*UBSdiff()/alphad_);

}


// ************************************************************************* //
