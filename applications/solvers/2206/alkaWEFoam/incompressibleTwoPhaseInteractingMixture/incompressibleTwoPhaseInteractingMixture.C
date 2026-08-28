/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2015 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

#include "incompressibleTwoPhaseInteractingMixture.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "fvc.H"
#include "../speciesProperties/speciesProperties.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(incompressibleTwoPhaseInteractingMixture, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct the gas-liquid mixture and the runtime-selected viscosity and
// capillary submodels used by the porous two-phase flow equations.
Foam::incompressibleTwoPhaseInteractingMixture::
incompressibleTwoPhaseInteractingMixture
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    const dictionary& dict
)
:
    twoPhaseMixture(U.mesh(), dict),
    dict_(dict),
    species1({"H2","O2","H2O"}),
	T_(
	        U.mesh().lookupObject<volScalarField>
        	(
            		"T"
        	)
	),
	MW_(
		U.mesh().lookupObject<speciesProperties>
        	(
            		"speciesProperties"
        	).MW()
	),

    
    
        muModel_
    (
        mixtureViscosityModel::New
        (
            "mu",
            dict.subDict(phase1Name_),
            dict.subDict(phase2Name_),
            U,
            phi
        )
    ),
    pCapModel_
    (
        capillaryPressureModel::New
        (
            alpha2_,
            dict_
        )
    ),
    nucModel_
    (
        viscosityModel::New
        (
            "nuc",
            dict.subDict(phase2Name_),
            U,
            phi
        )
    ),
    
        p_num_("p_num", dimPressure, dict),
    
    rhod_
    (
        IOobject
        (
            "rhod",
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        p_num_/(Foam::constant::physicoChemical::R*T_)*(MW_[2])
    ),
    rhoc_("rho", dimDensity, nucModel_->viscosityProperties()),
    
    dd_
    (
        "d",
        dimLength,
        muModel_->viscosityPropertiesSub1().getOrDefault("d", 0.0)
    ),
    alphaMax_(muModel_->viscosityPropertiesSub1().getOrDefault("alphaMax", 1.0)),

    U_(U),

    mu_
    (
        IOobject
        (
            "mu",
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U_.mesh(),
        dimensionedScalar(dimensionSet(1, -1, -1, 0, 0), Zero),
        calculatedFvPatchScalarField::typeName
    ),
    nu_
    (
        IOobject
        (
            "nu",
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mu_/rho()
    ),
  	   Pe_
       	(
    		U.mesh().lookupObject<volScalarField>
   	 	(
    			"Pe"
    			//(
    			//)
    		)
  	),
  	  	   Ne_
       	(
    		U.mesh().lookupObject<volScalarField>
   	 	(
    			"Ne"
    			//(
    			//)
    		)
  	),
  	  	   Mem_
       	(
    		U.mesh().lookupObject<volScalarField>
   	 	(
    			"Mem"
    			//(
    			//)
    		)
  	),
  	  	   PeC_
       	(
    		U.mesh().lookupObject<volScalarField>
   	 	(
    			"PeC"
    			//(
    			//)
    		)
  	),
  	  	   NeC_
       	(
    		U.mesh().lookupObject<volScalarField>
   	 	(
    			"NeC"
    			//(
    			//)
    		)
  	),
  	eps_
       	(
    		U.mesh().lookupObject<volScalarField>
   	 	(
    			"eps"
    			//(
    			//)
    		)
  	)
{

    
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::incompressibleTwoPhaseInteractingMixture::read()
{
    //if (regIOobject::read())
    //{
        if
        (
            muModel_().read(dict_.subDict(phase1Name_))
         && nucModel_().read(dict_.subDict(phase2Name_))
        )
        {
            muModel_->viscosityPropertiesSub1().readEntry("rho", rhod_);
            nucModel_->viscosityProperties().readEntry("rho", rhoc_);

            dd_ = dimensionedScalar
            (
                "d",
                dimLength,
                muModel_->viscosityPropertiesSub1().getOrDefault("d", 0)
            );

            alphaMax_ =
                muModel_->viscosityPropertiesSub1().getOrDefault
                (
                    "alphaMax",
                    1.0
                );

            return true;
        }
    //}

    return false;
}


// ************************************************************************* //
