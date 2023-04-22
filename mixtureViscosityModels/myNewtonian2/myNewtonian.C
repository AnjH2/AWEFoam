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

#include "myNewtonian.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace mixtureViscosityModels
{
    defineTypeNameAndDebug(myNewtonian, 0);

    addToRunTimeSelectionTable
    (
        mixtureViscosityModel,
        myNewtonian,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mixtureViscosityModels::myNewtonian::myNewtonian
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const word modelName
)
:
    mixtureViscosityModel(name, viscosityProperties, U, phi),
    myNewtonianCoeffs_(viscosityProperties.optionalSubDict(modelName + "Coeffs")),
    myNewtonianViscosityCoeff_("coeff", dimDynamicViscosity, myNewtonianCoeffs_),
    myNewtonianViscosityExponent_("exponent", dimless, myNewtonianCoeffs_),
    muMax_("muMax", dimDynamicViscosity, myNewtonianCoeffs_),
    alpha_
    (
        U.mesh().lookupObject<volScalarField>
        (
            IOobject::groupName
            (
                viscosityProperties.getOrDefault<word>("alpha", "alpha"),
                viscosityProperties.dictName()
            )
        )
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::mixtureViscosityModels::myNewtonian::mu(const volScalarField& muc) const
{
    return min
    (
        muc
      + myNewtonianViscosityCoeff_
       *(
            pow
            (
                scalar(10),
                myNewtonianViscosityExponent_*alpha_
            ) - scalar(1)
        ),
        muMax_
    );
}


bool Foam::mixtureViscosityModels::myNewtonian::read
(
    const dictionary& viscosityProperties
)
{
    mixtureViscosityModel::read(viscosityProperties);

    myNewtonianCoeffs_ = viscosityProperties.optionalSubDict(typeName + "Coeffs");

    myNewtonianCoeffs_.readEntry("k", myNewtonianViscosityCoeff_);
    myNewtonianCoeffs_.readEntry("n", myNewtonianViscosityExponent_);
    myNewtonianCoeffs_.readEntry("muMax", muMax_);

    return true;
}


// ************************************************************************* //
