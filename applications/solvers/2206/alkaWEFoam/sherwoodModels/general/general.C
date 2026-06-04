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

#include "general.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace sherwoodModels
{
    defineTypeNameAndDebug(general, 0);

    addToRunTimeSelectionTable
    (
        sherwoodModel,
        general,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sherwoodModels::general::general
(
    const word& name,
    const incompressibleTwoPhaseInteractingMixture& mixture,
    const dictionary& sherwoodPropertiesSub1,
    const volScalarField& lCh,
    const word modelName
)
:
    sherwoodModel(name, mixture, sherwoodPropertiesSub1,lCh),
    generalCoeffsSub1_(sherwoodPropertiesSub1.optionalSubDict(name+"_"+modelName + "Coeffs")),
    s_(generalCoeffsSub1_.getOrDefault<scalar>("scale",1)),
	abc_(3)
    

{
abc_.set
    	(
        	0,
        	new dimensionedScalar("a", dimless,generalCoeffsSub1_)
        );
abc_.set
    	(
        	1,
        	new dimensionedScalar("b", dimless,generalCoeffsSub1_)
        );
abc_.set
    	(
        	2,
        	new dimensionedScalar("c", dimless,generalCoeffsSub1_)
        );
        

        
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::sherwoodModels::general::ki(const int i)
	{
		return D2_[i]*sh(i)/(lCh_*s_);
	}



bool Foam::sherwoodModels::general::read
(
    const dictionary& sherwoodPropertiesSub1
)
{
    sherwoodModel::read(sherwoodPropertiesSub1);

    generalCoeffsSub1_ = sherwoodPropertiesSub1.optionalSubDict(typeName + "Coeffs");

    //generalCoeffsSub1_.readEntry("mu", mud_);
  //  generalCoeffs_.readEntry("n", generalsherwoodExponent_);
   // generalCoeffs_.readEntry("muMax", muMax_);

    return true;
}


// ************************************************************************* //
