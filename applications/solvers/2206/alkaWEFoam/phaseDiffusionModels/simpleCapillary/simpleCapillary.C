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

#include "simpleCapillary.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseDiffusionModels
{
    defineTypeNameAndDebug(simpleCapillary, 0);
    addToRunTimeSelectionTable(phaseDiffusionModel, simpleCapillary, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseDiffusionModels::simpleCapillary::simpleCapillary
(
    const dictionary& dict,
    const incompressibleTwoPhaseInteractingMixture& mixture,
    const word& modelName
)
:
    phaseDiffusionModel(dict, mixture,modelName),
    mixture_(mixture),
    eps_(
	    alphac_.mesh().lookupObject<volScalarField>
            (
            	"eps"
            )
        ),
    K_(
	    alphac_.mesh().lookupObject<volScalarField>
            (
            	"permeabilityField"
            )
        ),
    Solid_(
	    alphac_.mesh().lookupObject<volScalarField>
            (
            	"Solid"
            )
        )
{

}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseDiffusionModels::simpleCapillary::~simpleCapillary()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
volScalarField Foam::phaseDiffusionModels::simpleCapillary::A()
{

	return 1/K_*(mixture_.mud_m()/(krd()+SMALL)+mixture_.muc()/(krc()+SMALL));
}

volScalarField Foam::phaseDiffusionModels::simpleCapillary::krd()
{
	
	return pow(alphad_,mixture_.kr());
}
volScalarField Foam::phaseDiffusionModels::simpleCapillary::krc()
{
	
	return pow(1-alphad_,mixture_.kr());
}






// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::phaseDiffusionModels::simpleCapillary::correct()
{
    
    const volScalarField& dpcds(mixture_.dpcds());
    
 

    Ddm_=(-alphac_*rhoc_/mixture_.rho()) * 1/A()*dpcds*(1-Mem_)*tensor::I;
}


// ************************************************************************* //
