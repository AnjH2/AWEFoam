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

#include "simpleIshiiTEST.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace relativeVelocityModels
{
    defineTypeNameAndDebug(simpleIshiiTEST, 0);
    addToRunTimeSelectionTable(relativeVelocityModel, simpleIshiiTEST, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::relativeVelocityModels::simpleIshiiTEST::simpleIshiiTEST
(
    const dictionary& dict,
    const incompressibleTwoPhaseInteractingMixture& mixture
)
:
    relativeVelocityModel(dict, mixture),
    n_("n", dimless, dict),
    V0_("V0", dimVelocity, dict),
    V1_(
	alphac_.mesh().lookupObject<volScalarField>
        (
            	"V1"
        )
        ),
            eps_(
	alphac_.mesh().lookupObject<volScalarField>
        (
            	"eps"
        )
        ),
        D_("D", dimVelocity*dimLength, dict)

{

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::relativeVelocityModels::simpleIshiiTEST::~simpleIshiiTEST()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::relativeVelocityModels::simpleIshiiTEST::correct()
{

    Udm_ = (1-Mem_+VSMALL)*
    //(rhoc_/rho())*(V1_*V0_*pow(1-max(alphad_, scalar(0)),n_)+D_/alphad_*fvc::grad(alphad_));
    (rhoc_/rho())*(V1_*V0_*pow(1-max(alphad_, scalar(0)),n_));
    Udm_.correctBoundaryConditions();
    
    //Ddm_=(Pe_+Ne_)*(rhoc_/rho())*D_*mixture_.nu()/alphad_; old case <=1159
    //Ddm_=(Pe_+Ne_)*(rhoc_/rho())*D_/alphad_; 	case 1171 -> failed
    Ddm_=(1-Mem_+VSMALL)*(rhoc_/rho())*eps_*D_/alphad_;
    Ddm_.correctBoundaryConditions();
}


// ************************************************************************* //
