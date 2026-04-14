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

#include "sLinearVerticalVogt.H"
#include "addToRunTimeSelectionTable.H"
#include "../../porousProperties/porousProperties.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace coverageModels
{
    defineTypeNameAndDebug(sLinearVerticalVogt, 0);
    addToRunTimeSelectionTable(coverageModel, sLinearVerticalVogt, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coverageModels::sLinearVerticalVogt::sLinearVerticalVogt
(
    const dictionary& dict,
    const fvMesh& mesh,
    const incompressibleTwoPhaseInteractingMixture& mixture,
    const word& modelName 
)
:
    coverageModel(dict,mesh,mixture,modelName),
    electrodes({"Ne","Pe"}),
    dict_(dict),
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
		dimensionSet (0,-2,0,0,0,1,0),//find true dimension!!!--------------
		dict
	),
    J_scale_
	(
		"J_scale",
		dimless,//find true dimension!!!--------------
		dict
	),
	T_ref_
	(
		"Temperature_ref",
		dimTemperature,
		dict
	),
	CD_(electrodes.size()),
	c_(electrodes.size()),
	rd_(electrodes.size()),


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
    g_(meshObjects::gravity::New(mixture.U().time())),
    alphad_
    (
        mesh.lookupObject<volScalarField>
        (
            "alpha.gas"
        )
    )
    
{
forAll(electrodes,i)
	{
	CD_.set
    	(
        	i,
        	new dimensionedScalar("bubbleDragCoefficient_"+electrodes[i], dimless/(dimVelocity*dimVelocity),dict)
        );
    rd_.set
    	(
        	i,
        	new dimensionedScalar("rd_"+electrodes[i], dimLength,dict)
        );
     c_.set
        (
            i,
		    new dimensionedScalar("thetaScale_"+electrodes[i],dimless,dict)
	    );    
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::coverageModels::sLinearVerticalVogt::~sLinearVerticalVogt()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

volScalarField Foam::coverageModels::sLinearVerticalVogt::correction()
{

    return pow((1+(Pe_*CD_[1]+Ne_*CD_[0])*pow(mag(U_),2)),-2);
}
void Foam::coverageModels::sLinearVerticalVogt::correct()
{

    theta_ =max(alphad_/max((c_[0]*rd_[0]*as_[0]*Ne_+c_[1]*rd_[1]*as_[1]*Pe_),SMALL),correction()*(J_scale_*pow(mag(J_)/((Ne_*as_[0]+Pe_*as_[1]+as_[0]*VSMALL)*J_lim_),0.3))*pow(T_/T_ref_*dimensionedScalar(dimPressure,101325)/p_num_,2/3));//Numerical modeling and analysis of the effect of pressure on the performance of an alkaline water electrolysis system 

    theta_.correctBoundaryConditions();
}


// ************************************************************************* //
