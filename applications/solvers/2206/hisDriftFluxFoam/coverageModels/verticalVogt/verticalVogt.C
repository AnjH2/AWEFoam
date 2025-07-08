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

#include "verticalVogt.H"
#include "addToRunTimeSelectionTable.H"
#include "../../porousProperties/porousProperties.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace coverageModels
{
    defineTypeNameAndDebug(verticalVogt, 0);
    addToRunTimeSelectionTable(coverageModel, verticalVogt, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coverageModels::verticalVogt::verticalVogt
(
    const dictionary& dict,
    const fvMesh& mesh,
    const incompressibleTwoPhaseInteractingMixture& mixture,
    const word& modelName 
)
:
    coverageModel(dict,mesh,mixture,modelName),
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
	CD_
	(
		"bubbleDragCoeffient",
		dimless,
		dict
	),
	beta_
	(
		"bubbleContactAngle",
		dimless,
		dict
	),
	sigma_
	(
		"surfaceTension",
		dimForce/dimLength,
		dict
	),
	K2_
	(
		"correctionCoefficient",
		dimMass/dimLength/dimTime/dimTime,
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
    ),
    U_
    (
        mesh.lookupObject<volVectorField>
        (
            "U"
        )
    ),
    rhoc_(mixture_.rhoc()),
    rhod_(mixture_.rhod()),
    g_(meshObjects::gravity::New(mixture.U().time()))
    
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::coverageModels::verticalVogt::~verticalVogt()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

volScalarField Foam::coverageModels::verticalVogt::Cc()
{
    return  -6/(mag(g_)*(rhoc_-rhod_))*sin(degToRad()*beta_)*K2_;
}
volScalarField Foam::coverageModels::verticalVogt::Cb()
{
    return  3.0/4.0*CD_*rhoc_/(mag(g_)*(rhoc_-rhod_))*pow(mag(U_),2)*(1-(degToRad()*beta_-cos(degToRad()*beta_)*sin(degToRad()*beta_))/(Foam::constant::mathematical::pi));
}
volScalarField Foam::coverageModels::verticalVogt::correction()
{

    return pow((sqrt(-4*Cc()*dimensionedScalar(dimLength,1)+pow(Cb(),2))-Cb())/(2*sqrt(-Cc()*dimensionedScalar(dimLength,1))),4);
}
void Foam::coverageModels::verticalVogt::correct()
{
    //Info<<"Cb:"<<average(Cb())<<endl;
    
    theta_ =correction()*(J_scale_*pow(mag(J_)/((Ne_*as_[0]+Pe_*as_[1]+as_[0]*VSMALL)*J_lim_),0.3))*pow(T_/T_ref_*dimensionedScalar(dimPressure,101325)/p_num_,2/3);

    theta_.correctBoundaryConditions();
}


// ************************************************************************* //
