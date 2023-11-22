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

#include "WrongBeckermann.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace mixtureViscosityModels
{
    defineTypeNameAndDebug(WrongBeckermann, 0);

    addToRunTimeSelectionTable
    (
        mixtureViscosityModel,
        WrongBeckermann,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mixtureViscosityModels::WrongBeckermann::WrongBeckermann
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
    WrongBeckermannCoeffsSub1_(viscosityPropertiesSub1.optionalSubDict(modelName + "Coeffs")),
    WrongBeckermannCoeffsSub2_(viscosityPropertiesSub2.optionalSubDict(modelName + "Coeffs")),
    mud_("mu", dimDynamicViscosity, WrongBeckermannCoeffsSub1_),
    rhoc_("rho", dimDensity, WrongBeckermannCoeffsSub2_),
    //WrongBeckermannViscosityExponent_("exponent", dimless, WrongBeckermannCoeffs_),
    //muMax_("muMax", dimDynamicViscosity, WrongBeckermannCoeffs_),
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
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::mixtureViscosityModels::WrongBeckermann::mu(const volScalarField& muc, const volScalarField& rhod) const
{
    return pow(pow(1-alpha_,3)/(muc)+pow(alpha_,3)/(mud_),-1);
}


bool Foam::mixtureViscosityModels::WrongBeckermann::read
(
    const dictionary& viscosityPropertiesSub1
)
{
    mixtureViscosityModel::read(viscosityPropertiesSub1);

    WrongBeckermannCoeffsSub1_ = viscosityPropertiesSub1.optionalSubDict(typeName + "Coeffs");

    WrongBeckermannCoeffsSub1_.readEntry("mu", mud_);
  //  WrongBeckermannCoeffs_.readEntry("n", WrongBeckermannViscosityExponent_);
   // WrongBeckermannCoeffs_.readEntry("muMax", muMax_);

    return true;
}


// ************************************************************************* //
