/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2017 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "Beckermann.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace mixtureViscosityModels
{
    defineTypeNameAndDebug(Beckermann, 0);

    addToRunTimeSelectionTable
    (
        mixtureViscosityModel,
        Beckermann,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mixtureViscosityModels::Beckermann::Beckermann
(
    const word& name,
    const dictionary& viscosityPropertiesSub1,
    const dictionary& viscosityPropertiesSub2,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const word modelName
)
:
    mixtureViscosityModel(name, viscosityPropertiesSub1,viscosityPropertiesSub2, U, phi),
    BeckermannCoeffsSub1_(viscosityPropertiesSub1.optionalSubDict(modelName + "Coeffs")),
    BeckermannCoeffsSub2_(viscosityPropertiesSub2.optionalSubDict(modelName + "Coeffs")),
    mud_("mu", dimDynamicViscosity, BeckermannCoeffsSub1_),
    //BeckermannViscosityExponent_("exponent", dimless, BeckermannCoeffs_),
    //muMax_("muMax", dimDynamicViscosity, BeckermannCoeffs_),
    alpha_
    (
        U.mesh().lookupObject<volScalarField>
        (
            IOobject::groupName
            (
                viscosityPropertiesSub1.getOrDefault<word>("alpha", "alpha"),
                viscosityPropertiesSub1.dictName()
            )
        )
    ),
    rhod_
    (
        U.mesh().lookupObject<volScalarField>
        (
            IOobject::groupName
            (
                viscosityPropertiesSub1.getOrDefault<word>("rhod", "rhod"),
                viscosityPropertiesSub1.dictName()
            )
        )
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::mixtureViscosityModels::Beckermann::mu(const volScalarField& muc) const
{
    return pow(pow(alpha_,3)/muc+pow(1-alpha_,3)/mud_,-1);
}


bool Foam::mixtureViscosityModels::Beckermann::read
(
    const dictionary& viscosityPropertiesSub1
)
{
    mixtureViscosityModel::read(viscosityPropertiesSub1);

    BeckermannCoeffsSub1_ = viscosityPropertiesSub1.optionalSubDict(typeName + "Coeffs");

    BeckermannCoeffsSub1_.readEntry("mu", mud_);
  //  BeckermannCoeffs_.readEntry("n", BeckermannViscosityExponent_);
   // BeckermannCoeffs_.readEntry("muMax", muMax_);

    return true;
}


// ************************************************************************* //
