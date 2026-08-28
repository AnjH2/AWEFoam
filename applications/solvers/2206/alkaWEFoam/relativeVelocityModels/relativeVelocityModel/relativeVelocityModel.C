/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2017 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "relativeVelocityModel.H"
#include "fixedValueFvPatchFields.H"
#include "slipFvPatchFields.H"
#include "partialSlipFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(relativeVelocityModel, 0);
    defineRunTimeSelectionTable(relativeVelocityModel, dictionary);
}

// * * * * * * * * * * * * * Private Member Functions   * * * * * * * * * * * //

// Match Udm boundary treatment to the mixture-velocity boundary type so
// derived drift closures can update only the internal model expression.
Foam::wordList Foam::relativeVelocityModel::UdmPatchFieldTypes() const
{
    const volVectorField& U = mixture_.U();

    wordList UdmTypes
    (
        U.boundaryField().size(),
        calculatedFvPatchScalarField::typeName
    );

    forAll(U.boundaryField(), i)
    {
        if
        (
            isA<fixedValueFvPatchVectorField>(U.boundaryField()[i])
         || isA<slipFvPatchVectorField>(U.boundaryField()[i])
         || isA<partialSlipFvPatchVectorField>(U.boundaryField()[i])
        )
        {
            UdmTypes[i] = fixedValueFvPatchVectorField::typeName;
        }
    }

    return UdmTypes;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::relativeVelocityModel::relativeVelocityModel
(
    const dictionary& dict,
    const incompressibleTwoPhaseInteractingMixture& mixture,
    const word& modelName
)
:
    mixture_(mixture),
    dict_(dict),
    alphac_(mixture.alpha2()),
    alphad_(mixture.alpha1()),
    rhoc_(mixture.rhoc()),
    rhod_(mixture.rhod()),

    Udm_
    (
        IOobject
        (
            modelName+"Udm",
            alphac_.time().timeName(),
            alphac_.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        alphac_.mesh(),
        dimensionedVector(dimVelocity, Zero),
        UdmPatchFieldTypes()
    ),
    hF_(
    	dict.lookupOrDefault<scalar>("hF", 0) //how much the velocity is reduced, 0 is free rasing bubble, only active in Pe and Ne
    ),
        Pe_(
	        alphad_.mesh().lookupObject<volScalarField>
        	(
            		"Pe"
        	)
	),
	Ne_(
	        alphad_.mesh().lookupObject<volScalarField>
        	(
            		"Ne"
        	)
	),
	PeC_(
	        alphad_.mesh().lookupObject<volScalarField>
        	(
            		"PeC"
        	)
	),
	NeC_(
	        alphad_.mesh().lookupObject<volScalarField>
        	(
            		"NeC"
        	)
	),
	Mem_(
	        alphad_.mesh().lookupObject<volScalarField>
        	(
            		"Mem"
        	)
	),
		eps_(
	        alphad_.mesh().lookupObject<volScalarField>
        	(
            		"eps"
        	)
	),
	U_(
	        alphad_.mesh().lookupObject<volVectorField>
        	(
            		"U"
        	)
	),
	modelName_(modelName),
    // Phase dispersion is selected independently from the advective drift
    // closure, allowing buoyancy/lift and D.grad(alpha) to be treated separately.
    dModel_
    (
        phaseDiffusionModel::New
        (
            dict_,
            mixture,
            modelName_
        )
    )
{
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::relativeVelocityModel> Foam::relativeVelocityModel::New
(
    const dictionary& dict,
    const incompressibleTwoPhaseInteractingMixture& mixture,
    const word& modelName 
)
{
    //const word modelType(dict.get<word>(typeName));
    const word modelType(dict.get<word>(modelName+"RelativeVelocityModel"));
    Info<< "Selecting "<<modelName<<" relative velocity model " << modelType << endl;

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "relative velocity",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << abort(FatalIOError);
    }

    return
        autoPtr<relativeVelocityModel>
        (
            ctorPtr
            (
                dict.subDict(modelName+modelType + "Coeffs"),
                mixture,
                modelName
            )
        );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::relativeVelocityModel::~relativeVelocityModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::relativeVelocityModel::rho() const
{
    return alphac_*rhoc_ + alphad_*rhod_;
}


// Convert the gas relative velocity to the corresponding liquid velocity
// using the zero mass-weighted sum of phase-relative velocities.
Foam::tmp<Foam::volVectorField> Foam::relativeVelocityModel::Ucm() const
{
    volScalarField betac(alphac_*rhoc_);
    volScalarField betad(alphad_*rhod_);
    return -betad*Udm_/betac;
   
}


// Recover the drift velocity relative to the volume-averaged mixture from
// the mass-averaged phase-relative velocity stored in Udm_.
Foam::tmp<Foam::volVectorField> Foam::relativeVelocityModel::Udj() const
{
    return tmp<volVectorField>
    (
        new volVectorField
        (
            "Udj",
            rho()/rhoc_*Udm_
        )
    );
    
}

// Drift-induced momentum stress formed from the phase-relative velocities.
// Note: the present solver assembles its active drift stress directly from Ur.
Foam::tmp<Foam::volTensorField> Foam::relativeVelocityModel::tauDm() const
{

    return tmp<volTensorField>
    (
        new volTensorField
        (
            "tauDm",
            alphad_*rhod_*Udm_*Udm_+alphac_*rhoc_*Ucm()*Ucm()
        )
    );
}


// ************************************************************************* //
