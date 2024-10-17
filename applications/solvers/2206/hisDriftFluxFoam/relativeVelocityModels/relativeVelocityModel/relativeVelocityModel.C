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
    const incompressibleTwoPhaseInteractingMixture& mixture
)
:
    mixture_(mixture),
    alphac_(mixture.alpha2()),
    alphad_(mixture.alpha1()),
    rhoc_(mixture.rhoc()),
    rhod_(mixture.rhod()),

    Udm_
    (
        IOobject
        (
            "Udm",
            alphac_.time().timeName(),
            alphac_.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        alphac_.mesh(),
        dimensionedVector(dimVelocity, Zero),
        UdmPatchFieldTypes()
    ),
    Ddm_
    (
        IOobject
        (
            "Ddm",
            alphac_.time().timeName(),
            alphac_.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        alphac_.mesh(),
        dimensionedTensor(dimVelocity*dimLength, Zero)
    ),
    
    y0_(
    	dict.lookupOrDefault<scalar>("y0", Zero)//scalar(1.0)
    ),
    y1_(
    	dict.lookupOrDefault<scalar>("y1", Zero)
    ),
    y2_(
    	dict.lookupOrDefault<scalar>("y2", Zero)
    ),
    y3_(
    	dict.lookupOrDefault<scalar>("y3", Zero)
    ),
    yNormal_
    (
        IOobject
        (
            "yNormal",
            alphac_.time().timeName(),
            alphac_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alphac_.mesh(),
        dimensionedScalar(dimLength, 1)
    ),
    hF_(
    	dict.lookupOrDefault<scalar>("hF", 0) //how much the velocity is reduced, 0 is free rasing bubble, only active in Pe and Ne
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
	)
{
forAll ( alphac_.mesh().C(), celli) //loop through cell centres
{
  if(alphac_.mesh().C()[celli].y()<y1_) //not sure if this is correct syntax
  {
      if (mag(alphac_.mesh().C()[celli].y()-y0_)<=mag(alphac_.mesh().C()[celli].y()-y1_))
      {
      	yNormal_[celli]=alphac_.mesh().C()[celli].y()-y0_;
      }
      else
      {
      	yNormal_[celli]=alphac_.mesh().C()[celli].y()-y1_;
      }
  }
  else if ((alphac_.mesh().C()[celli].y()<y3_) and (alphac_.mesh().C()[celli].y()>y2_))
  {
      if (mag(alphac_.mesh().C()[celli].y()-y3_)<=mag(alphac_.mesh().C()[celli].y()-y2_))
      {
      	yNormal_[celli]=alphac_.mesh().C()[celli].y()-y3_;
      }
      else
      {
      	yNormal_[celli]=alphac_.mesh().C()[celli].y()-y2_;
      }
  }
  else
  {
  	yNormal_[celli]=1;
  }
}


}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::relativeVelocityModel> Foam::relativeVelocityModel::New
(
    const dictionary& dict,
    const incompressibleTwoPhaseInteractingMixture& mixture
)
{
    const word modelType(dict.get<word>(typeName));

    Info<< "Selecting relative velocity model " << modelType << endl;

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
                dict.optionalSubDict(modelType + "Coeffs"),
                mixture
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


// Calculate the relative velocity of the continuous phase w.r.t the mean
Foam::tmp<Foam::volVectorField> Foam::relativeVelocityModel::Ucm() const
{
    volScalarField betac(alphac_*rhoc_);
    volScalarField betad(alphad_*rhod_);
    return -betad*Udm_/betac;
    /*return tmp<volVectorField>
    (
        new volVectorField
        (
            "Ucm",
            -betad*Udm_/betac
        )
    );*/
    
}
// Calculate the relative velocity of the continuous phase w.r.t the mean
Foam::tmp<Foam::volTensorField> Foam::relativeVelocityModel::Dcm() const
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

// Calculate the relative velocity of the continuous phase w.r.t the mean
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
