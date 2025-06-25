/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2016 OpenFOAM Foundation
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

#include "vogt.H"
#include "addToRunTimeSelectionTable.H"
#include "../../porousProperties/porousProperties.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace coverageModels
{
    defineTypeNameAndDebug(vogt, 0);
    addToRunTimeSelectionTable(coverageModel, vogt, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coverageModels::vogt::vogt
(
    const dictionary& dict,
    const fvMesh& mesh,
    const word& modelName 
)
:
    coverageModel(dict,mesh,modelName),
    J_
    (
        mesh.lookupObject<volScalarField>
        (
            "J"
        )
    ),
    T_
    (
        mesh.lookupObject<volScalarField>
        (
            "T"
        )
    ),
        
    J_lim_
	(
		"J_lim",
		dimensionSet (0,-2,0,0,0,1,0),
		dict
	),
    J_scale_
	(
		"J_scale",
		dimless,
		dict
	),
	T_ref_
	(
		"Temperature_ref",
		dimTemperature,
		dict
	),
    as_(
	mesh.lookupObject<porousProperties>
        (
            "porousProperties"
        ).as()
        ),
    Pe_
    (
        mesh.lookupObject<volScalarField>
        (
            "Pe"
        )
    ),
    Ne_
    (
        mesh.lookupObject<volScalarField>
        (
            "Ne"
        )
    )
    
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::coverageModels::vogt::~vogt()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::coverageModels::vogt::correct()
{
    theta_ =(J_scale_*pow(mag(J_)/((Ne_*as_[0]+Pe_*as_[1]+as_[0]*VSMALL)*J_lim_),0.3))*pow(T_/T_ref_*dimensionedScalar(dimPressure,101325)/p_num_,2/3);//Numerical modeling and analysis of the effect of pressure on the performance of an alkaline water electrolysis system 

    theta_.correctBoundaryConditions();
}


// ************************************************************************* //
