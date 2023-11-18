/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2016 OpenFOAM Foundation
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

#include "fixedReferenceValue.H"
#include "addToRunTimeSelectionTable.H"
//#include "../../speciesMixture/speciesMixture.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace concentrationLossModels
{
    defineTypeNameAndDebug(fixedReferenceValue, 0);
    addToRunTimeSelectionTable(concentrationLossModel, fixedReferenceValue, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::concentrationLossModels::fixedReferenceValue::fixedReferenceValue
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    concentrationLossModel(dict,mesh),
    C2_Ref_(species2.size()),
    n_Ne_(species2.size()),
    n_Pe_(species2.size()),
    C2_s_H2_
   	(
    		mesh.lookupObject<volScalarField>
   	 	(
    			"C2_s_H2"
    			//(
    			//)
    		)
  	),
    C2_s_O2_
   	(
    		mesh.lookupObject<volScalarField>
   	 	(
    			"C2_s_O2"
    			//(
    			//)
    		)
  	),
    C2_s_H2O_
   	(
    		mesh.lookupObject<volScalarField>
   	 	(
    			"C2_s_H2O"
    			//(
    			//)
    		)
  	),
    C2_s_OH_
   	(
    		mesh.lookupObject<volScalarField>
   	 	(
    			"C2_s_OH"
    			//(
    			//)
    		)
  	)
    /*C2_H2_
   	(
    		mesh.lookupObject<volScalarField>
   	 	(
    			"C_H2.electrolyte"
    			//(
    			//)
    		)
  	),
        C2_O2_
   	(
    		mesh.lookupObject<volScalarField>
   	 	(
    			"C_O2.electrolyte"
    			//(
    			//)
    		)
  	),
  	    C2_H2O_
   	(
    		mesh.lookupObject<volScalarField>
   	 	(
    			"C_H2O.electrolyte"
    			//(
    			//)
    		)
  	),
  	    C2_OH_
   	(
    		mesh.lookupObject<volScalarField>
   	 	(
    			"C_OH.electrolyte"
    			//(
    			//)
    		)
  	)*/
    //C2_(species2.size())
    /*C2_(
		mesh.lookupObject<speciesMixture>
        	(
        	    "speciesMixture"
        	).C2()
        )*/
{
forAll(species2,i)
	{
	n_Ne_.set
    	(
        	i,
        	new dimensionedScalar("n_Ne_"+species2[i], dimless,dict)
        );
        n_Pe_.set
    	(
        	i,
        	new dimensionedScalar("n_Pe_"+species2[i], dimless,dict)
        );
	C2_Ref_.set
        	(
        	i,
        	new dimensionedScalar("C_ref_"+species2[i], dimensionSet(0, -3, 0, 0, 1, 0, 0),dict)
        	);   	
	}
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::concentrationLossModels::fixedReferenceValue::~fixedReferenceValue()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::concentrationLossModels::fixedReferenceValue::correct(const PtrList <volScalarField>& C2_s_)
{
	forAll(species2,i)
	{
	CR_Ne_[i]=Foam::pow(C2_s_[i]/C2_Ref_[i],n_Ne_[i]);
	CR_Ne_[i].correctBoundaryConditions();
	CR_Pe_[i]=Foam::pow(C2_s_[i]/C2_Ref_[i],n_Pe_[i]);
	CR_Pe_[i].correctBoundaryConditions();
	}
}


// ************************************************************************* //
