/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "speciesMixture.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::speciesMixture::speciesMixture
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phase1NamE_,
    const word& phase2NamE_
)
:
    //phase1NamE_(wordList(dict.lookup("phases"))[0]),
    //phase2NamE_(wordList(dict.lookup("phases"))[1]),

    C_H2_1_
    (
        IOobject
        (
            IOobject::groupName("C_H2_", phase1NamE_),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    
    C_O2_1_
    (
        IOobject
        (
            IOobject::groupName("C_O2_", phase1NamE_),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    
    C_OH_1_
    (
        IOobject
        (
            IOobject::groupName("C_OH_", phase1NamE_),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    
    C_H2O_1_
    (
        IOobject
        (
            IOobject::groupName("C_H2O_", phase1NamE_),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    C_H2_2_
    (
        IOobject
        (
            IOobject::groupName("C_H2_", phase2NamE_),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    
    C_O2_2_
    (
        IOobject
        (
            IOobject::groupName("C_O2_", phase2NamE_),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    )
{}


// ************************************************************************* //
