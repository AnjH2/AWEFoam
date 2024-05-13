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

#include "Whitaker.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace sherwoodModels
{
    defineTypeNameAndDebug(Whitaker, 0);

    addToRunTimeSelectionTable
    (
        sherwoodModel,
        Whitaker,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sherwoodModels::Whitaker::Whitaker
(
    const word& name,
    const incompressibleTwoPhaseInteractingMixture& mixture,
    const dictionary& sherwoodPropertiesSub1,
    const word modelName
)
:
    sherwoodModel(name, mixture, sherwoodPropertiesSub1),
    WhitakerCoeffsSub1_(sherwoodPropertiesSub1.optionalSubDict(name+"_"+modelName + "Coeffs")),
        d_
	(
		"d",
		dimensionSet ( 0, 1, 0, 0, 0, 0,0),
		WhitakerCoeffsSub1_
	)
    

{
Info<<"\t" <<modelName << " model should only be used for moving bubbles"<<endl;
Info<<"\tAssuming dynamic viscosity dosent change significantly as function of concentration"<<endl;

}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::sherwoodModels::Whitaker::ki(const int i)
	{
		return D2_[i]*sh(i)/(d_);
	}



bool Foam::sherwoodModels::Whitaker::read
(
    const dictionary& sherwoodPropertiesSub1
)
{
    sherwoodModel::read(sherwoodPropertiesSub1);

    WhitakerCoeffsSub1_ = sherwoodPropertiesSub1.optionalSubDict(typeName + "Coeffs");

    //WhitakerCoeffsSub1_.readEntry("mu", mud_);
  //  WhitakerCoeffs_.readEntry("n", WhitakersherwoodExponent_);
   // WhitakerCoeffs_.readEntry("muMax", muMax_);

    return true;
}


// ************************************************************************* //
