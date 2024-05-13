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


#include "sherwoodModel.H"
#include "volFields.H"
#include "surfaceMesh.H"
#include "../../speciesProperties/speciesProperties.H"
#include "../../relativeVelocityModels/relativeVelocityModel/relativeVelocityModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sherwoodModel, 0);
    defineRunTimeSelectionTable(sherwoodModel, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sherwoodModel::sherwoodModel
(
    const word& name,
    const incompressibleTwoPhaseInteractingMixture& mixture,
    const dictionary& sherwoodPropertiesSub1
)
:
    name_(name),
    mixture_(mixture),
    species2({"H2","O2","H2O","OH"}),
    sherwoodPropertiesSub1_(sherwoodPropertiesSub1),

    D2_(
    	mixture_.rhod().mesh().lookupObject<speciesProperties>
        	(
            		"speciesProperties"
        	).D2()
        ),
            Udm_(
    	mixture_.rhod().mesh().lookupObject<volVectorField>
        	(
            		"Udm"
        	)
        )
	
{

}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
bool Foam::sherwoodModel::read(const dictionary& sherwoodPropertiesSub1)
{
    sherwoodPropertiesSub1_ = sherwoodPropertiesSub1;

    return true;
}

// ************************************************************************* //
