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
    along with OpenFOAM.  If not, see <http://(www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "speciesTransport.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "../reactionProperties/reactionProperties.H"
#include "../porousProperties/porousProperties.H"
#include "../speciesProperties/speciesProperties.H"
#include "../incompressibleTwoPhaseInteractingMixture/incompressibleTwoPhaseInteractingMixture.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::speciesTransport::speciesTransport
(
    const fvMesh& mesh,
    const dictionary& dict,
    const massAndSpeciesTransferModel& mSTaPtr
)
:
	mSTaPtr_(mSTaPtr),
    	speciesTransportCoeffs_(dict.optionalSubDict("speciesTransport")),
    shModelW_
    (
        sherwoodModel::New
        (
            "wireModel",
            mSTaPtr_.mixture(),
            speciesTransportCoeffs_
        )
    ),
    	species2({"H2","O2","H2O","OH"}),
    	C2_s_(species2.size()),
	C2_(
		mesh.lookupObject<speciesProperties>
        	(
            		"speciesProperties"
        	).C2()
	),
	alphal_(
	        mesh.lookupObject<volScalarField>
        	(
            		"alpha.electrolyte"
        	)
	),
		Pe_(
	        mesh.lookupObject<volScalarField>
        	(
            		"Pe"
        	)
	),
	Ne_(
	        mesh.lookupObject<volScalarField>
        	(
            		"Ne"
        	)
	),
	Mem_(
	        mesh.lookupObject<volScalarField>
        	(
            		"Mem"
        	)
	),

	as_(
		mesh.lookupObject<porousProperties>
        	(
            		"porousProperties"
        	).as()
	),

    	MW_(
		mesh.lookupObject<speciesProperties>
        	(
            		"speciesProperties"
        	).MW()
        )

{
forAll(species2,i)
	{
        	C2_s_.set
        	(
        	        i,
        		new volScalarField
        		(
            			IOobject
            			(
               				"C2_s_"+species2[i],
                			mesh.time().timeName(),
                			mesh,
                			IOobject::NO_READ,
                			IOobject::NO_WRITE
            			),
            			C2_[i]
        		)
        	);
	}
}
// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::speciesTransport::correct(const int i, const volScalarField& theta_)
{

		C2_s_[i]=(mSTaPtr_.Psi_BV()[i])/(shModelW_->ki(i)*as_[0]*(1-theta_))+C2_[i]; //original
		//C2_s_[i]=(mSTaPtr_.Psi_BV()[i])/(shModelW_->ki(i)*as_[0]*(1-theta_))+C2_[i]; // does not have the error
		//C2_s_[i]=(mSTaPtr_.Psi_BV()[i]-mSTaPtr_.mDot_Wall()[i]/MW_[i])/(shModelW_->ki(i)*as_[0])+C2_[i]; // have less error
}



// ************************************************************************* //
