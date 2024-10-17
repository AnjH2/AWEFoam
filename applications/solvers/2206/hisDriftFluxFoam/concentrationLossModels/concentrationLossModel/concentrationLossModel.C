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
    const fvMesh& mesh,
    const massAndSpeciesTransferModel& mSTaPtr
)
:
	concentrationLossModelDict_(dict),
	species2({"H2","O2","H2O","OH"}),
	CR_(species2.size()),
    	CRa_Ne_(
    		IOobject
    		(
    	    		"CRa_Ne",
            		mesh.time().timeName(),
            		mesh,
            		IOobject::NO_READ,
            		IOobject::AUTO_WRITE
    		),
    		mesh,
    		dimensionedScalar(dimless,One)
    	),
    	CRa_Pe_(
    		IOobject
    		(
    	    		"CRa_Pe",
            		mesh.time().timeName(),
            		mesh,
            		IOobject::NO_READ,
            		IOobject::AUTO_WRITE
    		),
    		mesh,
    		dimensionedScalar(dimless,One)
    	),
    	CRc_Ne_(
    		IOobject
    		(
    	    		"CRc_Ne",
            		mesh.time().timeName(),
            		mesh,
            		IOobject::NO_READ,
            		IOobject::AUTO_WRITE
    		),
    		mesh,
    		dimensionedScalar(dimless,One)
    	),
    	CRc_Pe_(
    		IOobject
    		(
    	    		"CRc_Pe",
            		mesh.time().timeName(),
            		mesh,
            		IOobject::NO_READ,
            		IOobject::AUTO_WRITE
    		),
    		mesh,
    		dimensionedScalar(dimless,One)
    	),
    	Pe_
    	(
    		mesh.lookupObject<volScalarField>
        	(
            		"Pe"
        	)
        ),
    	Ne_
    	(
    		mesh.lookupObject<volScalarField>
        	(
            		"Ne"
        	)
       	),
    	a_Ne_(concentrationLossModelDict_.get<wordList>("a_Ne")),
	a_Pe_(concentrationLossModelDict_.get<wordList>("a_Pe")),
	c_Ne_(concentrationLossModelDict_.get<wordList>("c_Ne")),
	c_Pe_(concentrationLossModelDict_.get<wordList>("c_Pe")),
	ab_Ne_(species2.size()),
	ab_Pe_(species2.size()),
	cb_Ne_(species2.size()),
	cb_Pe_(species2.size())/*,
	VOne_(
    		IOobject
    		(
    	    		"VOne",
            		mesh.time().timeName(),
            		mesh,
            		IOobject::NO_READ,
            		IOobject::AUTO_WRITE
    		),
    		mesh,
    		dimensionedScalar(dimless,One)
    	)*/
{
forAll(species2,i)
	{
	forAll(a_Ne_,j){
		if (a_Ne_[j]==species2[i]){
			ab_Ne_[i]=1;
			break;
		} else {
			ab_Ne_[i]=0;
		}	
	}
		forAll(a_Pe_,j){
		if (a_Pe_[j]==species2[i]){
			ab_Pe_[i]=1;
			break;
		} else {
			ab_Pe_[i]=0;
		}
		
	}
		forAll(c_Ne_,j){
		if (c_Ne_[j]==species2[i]){
			cb_Ne_[i]=1;
			break;
		} else {
			cb_Ne_[i]=0;
		}
		
	}
		forAll(c_Pe_,j){
		if (c_Pe_[j]==species2[i]){
			cb_Pe_[i]=1;
			break;
		} else {
			cb_Pe_[i]=0;
		}
		
	}

	CR_.set
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
            		dimensionedScalar("CR_"+species2[i],dimless,One)
        	)
        );
	}
	/*
	Info<<"ab_Ne"<<ab_Ne_<<endl;
	Info<<"ab_Pe"<<ab_Pe_<<endl;
	Info<<"cb_Ne"<<cb_Ne_<<endl;
	Info<<"cb_Pe"<<cb_Pe_<<endl;
	*/
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::concentrationLossModel> Foam::concentrationLossModel::New
(
    const dictionary& dict,
    const fvMesh& mesh,
    const massAndSpeciesTransferModel& mSTaPtr
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
                mesh,
                mSTaPtr
            )
        );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::concentrationLossModel::~concentrationLossModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //




// ************************************************************************* //
