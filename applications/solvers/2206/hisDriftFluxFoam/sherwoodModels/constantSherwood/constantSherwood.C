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

#include "constantSherwood.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace sherwoodModels
{
    defineTypeNameAndDebug(constantSherwood, 0);

    addToRunTimeSelectionTable
    (
        sherwoodModel,
        constantSherwood,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sherwoodModels::constantSherwood::constantSherwood
(
    const word& name,
    const incompressibleTwoPhaseInteractingMixture& mixture,
    const dictionary& sherwoodPropertiesSub1,
    const word modelName
)
:
    sherwoodModel(name, mixture, sherwoodPropertiesSub1),
    constantSherwoodCoeffsSub1_(sherwoodPropertiesSub1.optionalSubDict(name+"_"+modelName + "Coeffs")),
        d_
	(
		"d",
		dimensionSet ( 0, 1, 0, 0, 0, 0,0),
		constantSherwoodCoeffsSub1_
	),
    sh_(species2.size())
    

{
forAll(species2,i)
	{
        sh_.set
    	(
        	i,
        	new dimensionedScalar("sh_"+species2[i], dimless,constantSherwoodCoeffsSub1_)
        );

        
	}

}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::sherwoodModels::constantSherwood::ki(const int i)
	{
		return D2_[i]*sh_[i]/(d_);
	}

bool Foam::sherwoodModels::constantSherwood::read
(
    const dictionary& sherwoodPropertiesSub1
)
{
    sherwoodModel::read(sherwoodPropertiesSub1);

    constantSherwoodCoeffsSub1_ = sherwoodPropertiesSub1.optionalSubDict(typeName + "Coeffs");

    //constantSherwoodCoeffsSub1_.readEntry("mu", mud_);
  //  constantSherwoodCoeffs_.readEntry("n", constantSherwoodsherwoodExponent_);
   // constantSherwoodCoeffs_.readEntry("muMax", muMax_);

    return true;
}


// ************************************************************************* //
