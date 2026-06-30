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
    under the terms of the GNU generalR Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU generalR Public License
    for more details.

    You should have received a copy of the GNU generalR Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "generalR.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace sherwoodModels
{
    defineTypeNameAndDebug(generalR, 0);

    addToRunTimeSelectionTable
    (
        sherwoodModel,
        generalR,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sherwoodModels::generalR::generalR
(
    const word& name,
    const incompressibleTwoPhaseInteractingMixture& mixture,
    const dictionary& sherwoodPropertiesSub1,
    const word modelName
)
:
    sherwoodModel(name, mixture, sherwoodPropertiesSub1),
    generalRCoeffsSub1_(sherwoodPropertiesSub1.optionalSubDict(name+"_"+modelName + "Coeffs")),
        d_(electrodes.size()),
	abc_(3)
    

{
abc_.set
    	(
        	0,
        	new dimensionedScalar("a", dimless,generalRCoeffsSub1_)
        );
abc_.set
    	(
        	1,
        	new dimensionedScalar("b", dimless,generalRCoeffsSub1_)
        );
abc_.set
    	(
        	2,
        	new dimensionedScalar("c", dimless,generalRCoeffsSub1_)
        );
forAll(electrodes,i)
	{
	d_.set
    	(
        	i,
        	new dimensionedScalar("d_"+electrodes[i], dimLength,generalRCoeffsSub1_)
        );

        }
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::sherwoodModels::generalR::ki(const int i)
	{
		return D2_[i]*sh(i)/(d_[0]*(Ne_+NeC_)+d_[1]*(Pe_+PeC_)+(d_[1]+d_[0])/2*Mem_);
	}



bool Foam::sherwoodModels::generalR::read
(
    const dictionary& sherwoodPropertiesSub1
)
{
    sherwoodModel::read(sherwoodPropertiesSub1);

    generalRCoeffsSub1_ = sherwoodPropertiesSub1.optionalSubDict(typeName + "Coeffs");

    //generalRCoeffsSub1_.readEntry("mu", mud_);
  //  generalRCoeffs_.readEntry("n", generalRsherwoodExponent_);
   // generalRCoeffs_.readEntry("muMax", muMax_);

    return true;
}


// ************************************************************************* //
