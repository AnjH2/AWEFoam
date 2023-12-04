/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2017 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "concentrationLossModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(concentrationLossModel, 0);
    defineRunTimeSelectionTable(concentrationLossModel, dictionary);
}

// * * * * * * * * * * * * * Private Member Functions   * * * * * * * * * * * //




// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::concentrationLossModel::concentrationLossModel
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
	concentrationLossModelDict_(dict),
	species2({"H2","O2","H2O","OH"}),
    	CR_Ne_(species2.size()),
    	CR_Pe_(species2.size())

{
forAll(species2,i)
	{
	CR_Ne_.set
    	(
        	i,
        	new volScalarField
        	(
        		IOobject
        		(
        			"CR_"+species2[i],
        			mesh.time().timeName(),
                		mesh,
                		IOobject::NO_READ,
                		IOobject::AUTO_WRITE
            		),
            		mesh,
            		dimensionedScalar("CR__"+species2[i],dimless,One)
        	)
        );
        CR_Pe_.set
    	(
        	i,
        	new volScalarField
        	(
        		IOobject
        		(
        			"CR_"+species2[i],
        			mesh.time().timeName(),
                		mesh,
                		IOobject::NO_READ,
                		IOobject::AUTO_WRITE
            		),
            		mesh,
            		dimensionedScalar("CR__"+species2[i],dimless,One)
        	)
        );	
	}
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::concentrationLossModel> Foam::concentrationLossModel::New
(
    const dictionary& dict,
    const fvMesh& mesh
)
{
    const word modelType(dict.get<word>(typeName));

    Info<< "Selecting concentrationLoss model " << modelType << endl;

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "concentrationLoss model",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << abort(FatalIOError);
    }

    return
        autoPtr<concentrationLossModel>
        (
            ctorPtr
            (
                dict.optionalSubDict(modelType + "Coeffs"),
                mesh
            )
        );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::concentrationLossModel::~concentrationLossModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //




// ************************************************************************* //
