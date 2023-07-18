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

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(incompressibleTwoPhaseInteractingMixture, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::incompressibleTwoPhaseInteractingMixture::
incompressibleTwoPhaseInteractingMixture
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    IOdictionary
    (
        IOobject
        (
            "transportProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    twoPhaseMixture(U.mesh(), *this),
    speciesMixture(U.mesh(), *this, phase1Name_, phase2Name_),
    
    muModel_
    (
        mixtureViscosityModel::New
        (
            "mu",
            subDict(phase1Name_),
            U,
            phi
        )
    ),

    nucModel_
    (
        viscosityModel::New
        (
            "nuc",
            subDict(phase2Name_),
            U,
            phi
        )
    ),

    //rhod_("rho", dimDensity, muModel_->viscosityProperties()),
    
    
    p_num_("p_num", dimPressure, *this),
    
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
        p_num_/(Foam::constant::physicoChemical::R*T_)*(C1_[0]*MW_[0]+C1_[1]*MW_[1])*(1/(C1_[0]+C1_[1]))
    ),
    rhoc_("rho", dimDensity, nucModel_->viscosityProperties()),
    
    /*rhoc_
    (
        IOobject
        (
            "rhoc",
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        C2_[0]*MW_[0]+C2_[1]*MW_[1]+C2_[2]*MW_[2]+C2_[3]*MW_[3]+C2_[3]*MW_[4]
    ),*/
    dd_
    (
        "d",
        dimLength,
        muModel_->viscosityProperties().getOrDefault("d", 0.0)
    ),
    alphaMax_(muModel_->viscosityProperties().getOrDefault("alphaMax", 1.0)),

    U_(U),
    phi_(phi),

    mu_
    (
        IOobject
        (
            "mu",
            U_.time().timeName(),
            U_.db()
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
            IOobject::NO_WRITE
        ),
        mu_/rho()
    ),
    
    prgh_
       	(
    		U.mesh().lookupObject<volScalarField>
   	 	(
    			"p_rgh"
    			//(
    			//)
    		)
  	)
{
    correct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::incompressibleTwoPhaseInteractingMixture::read()
{
    if (regIOobject::read())
    {
        if
        (
            muModel_().read(subDict(phase1Name_))
         && nucModel_().read(subDict(phase2Name_))
        )
        {
            muModel_->viscosityProperties().readEntry("rho", rhod_);
            nucModel_->viscosityProperties().readEntry("rho", rhoc_);

            dd_ = dimensionedScalar
            (
                "d",
                dimLength,
                muModel_->viscosityProperties().getOrDefault("d", 0)
            );

            alphaMax_ =
                muModel_->viscosityProperties().getOrDefault
                (
                    "alphaMax",
                    1.0
                );

            return true;
        }
    }

    return false;
}
/*void Foam::incompressibleTwoPhaseInteractingMixture::correct()
{
	mu_ = muModel_->mu(rhoc_*nucModel_->nu());
	D1_eff_=D1_ref_*(1+DeltaX1(C1_[1]/(C1_[0]+C1_[1])))/(1+DeltaX1(dummy_));
}*/

// ************************************************************************* //
