/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014 OpenFOAM Foundation
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

#include "capillaryPressureModel.H"
#include "volFields.H"
#include "surfaceMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(capillaryPressureModel, 0);
    defineRunTimeSelectionTable(capillaryPressureModel, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::capillaryPressureModel::capillaryPressureModel
(
    const volScalarField& alphaWetting,
    const dictionary& dict
)
:
    alphaWetting_(alphaWetting),
    dict_(dict),
    Solid_
    (
    		alphaWetting_.mesh().lookupObject<volScalarField>
   	 	    (
    			"Solid"
    			//(
    			//)
    		)
  	),
  	pc_
    (
        IOobject
        (
            "pc",
            alphaWetting_.time().timeName(),
            alphaWetting_.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),       
        alphaWetting_.mesh(),
        dimensionedScalar("pc", dimPressure, 0)
    ),
    dpcds_
    (
        IOobject
        (
            "dpcds",
            alphaWetting_.time().timeName(),
            alphaWetting_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),       
        alphaWetting_.mesh(),
        dimensionedScalar("dpcds", dimPressure, 0)
    )
   
{

}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::capillaryPressureModel::read(const dictionary& dict)
{
    const dictionary& dict_ = dict;

    return true;
}


// ************************************************************************* //
