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

#include "phaseDiffusionModel.H"
#include "fixedValueFvPatchFields.H"
#include "slipFvPatchFields.H"
#include "partialSlipFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(phaseDiffusionModel, 0);
    defineRunTimeSelectionTable(phaseDiffusionModel, dictionary);
}

// * * * * * * * * * * * * * Private Member Functions   * * * * * * * * * * * //

Foam::wordList Foam::phaseDiffusionModel::UdmPatchFieldTypes() const
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

Foam::phaseDiffusionModel::phaseDiffusionModel
(
    const dictionary& dict,
    const incompressibleTwoPhaseInteractingMixture& mixture,
    const word& modelName
)
:
    mixture_(mixture),
    alphac_(mixture.alpha2()),
    alphad_(mixture.alpha1()),
    rhoc_(mixture.rhoc()),
    rhod_(mixture.rhod()),
    Ddm_
    (
        IOobject
        (
            modelName+"Ddm",
            alphac_.time().timeName(),
            alphac_.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        alphac_.mesh(),
        dimensionedTensor(dimVelocity*dimLength, Zero)
    ),
    dF_(
    	dict.lookupOrDefault<scalar>("dF", 1) //how much the dispersion is incressed, 1 is defined dispersion, only active in Pe and Ne
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
	modelName_(modelName)
{

}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::phaseDiffusionModel> Foam::phaseDiffusionModel::New
(
    const dictionary& dict,
    const incompressibleTwoPhaseInteractingMixture& mixture,
    const word& modelName 
)
{
    //const word modelType(dict.get<word>(typeName));
    const word modelType(dict.get<word>(modelName+"phaseDiffusionModel"));
    Info<< "Selecting phase diffusion model " << modelType << endl;

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "phase diffusion model",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << abort(FatalIOError);
    }

    return
        autoPtr<phaseDiffusionModel>
        (
            ctorPtr
            (
                dict.subDict(modelType + "Coeffs"),
                mixture,
                modelName
            )
        );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseDiffusionModel::~phaseDiffusionModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
volScalarField Foam::phaseDiffusionModel::VToM()
{

	return rhoc_/mixture_.rho();//from gas volume center to gas mass center
}
// Calculate the relative velocity of the continuous phase w.r.t the mean
Foam::tmp<Foam::volTensorField> Foam::phaseDiffusionModel::Dcm() const
{
    volScalarField betac(alphac_*rhoc_);
    volScalarField betad(alphad_*rhod_);
    return tmp<volTensorField>
    (
        new volTensorField
        (
            "Dcm",
            -betad*Ddm_/betac
        )
    );
    
}



// ************************************************************************* //
