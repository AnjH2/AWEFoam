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

#include "localReference.H"
#include "addToRunTimeSelectionTable.H"
#include "../../reactionProperties/reactionProperties.H"
#include "../../porousProperties/porousProperties.H"
#include "../../speciesProperties/speciesProperties.H"
#include "../../speciesTransport/speciesTransport.H"
#include "../../incompressibleTwoPhaseInteractingMixture/incompressibleTwoPhaseInteractingMixture.H"
//#include "../../speciesMixture/speciesMixture.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace concentrationLossModels
{
    defineTypeNameAndDebug(localReference, 0);
    addToRunTimeSelectionTable(concentrationLossModel, localReference, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::concentrationLossModels::localReference::localReference
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    concentrationLossModel(dict,mesh),
    n_Ne_(species2.size()),
    n_Pe_(species2.size()),
    /*C2_s_(
	mesh.lookupObject<speciesTransport>
        (
            "transportProperties"
        ).C2_s()
        ),*/
    C2_(
	mesh.lookupObject<speciesProperties>
        (
            "speciesProperties"
        ).C2()
        ),
    nb_(concentrationLossModelDict_.get<bool>("stoichiometricCoefficient"))
    
{
if (nb_) {
	forAll(species2,i) {
		n_Ne_[i]=
		(
			mesh.lookupObject<reactionProperties>
        		(
            			"reactionProperties"
        		).psi()[i]
        	);
        	n_Pe_[i]=
		(
			mesh.lookupObject<reactionProperties>
        		(
            			"reactionProperties"
        		).psi()[i+4]
        	);
        }
} else {
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
        }
}
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::concentrationLossModels::localReference::~localReference()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::concentrationLossModels::localReference::correct(const PtrList <volScalarField>& C2_s_)
{
	
	CRa_Ne_=VOne_;
	CRa_Pe_=VOne_;
	CRc_Ne_=VOne_;
	CRc_Pe_=VOne_;
	forAll(species2,i)
	{
		CR_[i]=C2_s_[i]/C2_[i];
		if (ab_Ne_[i]==1) {
			CRa_Ne_=CRa_Ne_*Foam::pow(CR_[i],n_Ne_[i]);
		}
		if (ab_Pe_[i]==1) {
			CRa_Pe_=CRa_Pe_*Foam::pow(CR_[i],n_Pe_[i]);
		}
		if (cb_Ne_[i]==1) {
			CRc_Ne_=CRc_Ne_*Foam::pow(CR_[i],n_Ne_[i]);
		}
		if (cb_Pe_[i]==1) {
			CRc_Pe_=CRc_Pe_*Foam::pow(CR_[i],n_Pe_[i]);
		}
	}
}


// ************************************************************************* //
