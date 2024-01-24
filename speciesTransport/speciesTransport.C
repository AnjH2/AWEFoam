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
    	species2({"H2","O2","H2O","OH"}),
    	C2_s_(species2.size()),
    	k_as_(species2.size()),
    	Re_p_(
            	IOobject
            	(
               		"Re_p_",
                	mesh.time().timeName(),
                	mesh,
                	IOobject::NO_READ,
                	IOobject::NO_WRITE
            	),
            	mesh,
            	dimless
        
        ),
        ReSMALLF_(
            	IOobject
            	(
               		"ReSMALLF_",
                	mesh.time().timeName(),
                	mesh,
                	IOobject::NO_READ,
                	IOobject::NO_WRITE
            	),
            	mesh,
            	dimensionedScalar(dimless,1)
        
        ),
    	Sc_(species2.size()),
    	SH_(species2.size()),
	sh_(3),
	D2_(
		mesh.lookupObject<speciesProperties>
        	(
            		"speciesProperties"
        	).D2()
	),
	C2_(
		mesh.lookupObject<speciesProperties>
        	(
            		"speciesProperties"
        	).C2()
	),
	/*Psi_BV_(
		mesh.lookupObject<massAndSpeciesTransferModel>
        	(
            		"transportProperties"
        	).Psi_BV()
	),*/
	/*Psi_(
		mesh.lookupObject<massAndSpeciesTransferModel>
        	(
            		"transportProperties"
        	).Psi()
	
	),*/
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
	nuc_(
	        mesh.lookupObject<volScalarField>
        	(
            		"nuc"
        	)
	),
	U_(
	        mesh.lookupObject<volVectorField>
        	(
            		"U"
        	)
	),
	as_(
		mesh.lookupObject<porousProperties>
        	(
            		"porousProperties"
        	).as()
	),
	D_pore_(
		mesh.lookupObject<porousProperties>
        	(
            		"porousProperties"
        	).D_pore()
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
                			IOobject::AUTO_WRITE
            			),
            			mesh,
            			dimensionSet(0,-3,0,0,1,0,0)
        		)
        	);
        	k_as_.set
        	(
        	        i,
        		new volScalarField
        		(
            			IOobject
            			(
               				"k_as_"+species2[i],
                			mesh.time().timeName(),
                			mesh,
                			IOobject::NO_READ,
                			IOobject::NO_WRITE
            			),
            			mesh,
            			dimensionSet(0,1,-1,0,0,0,0)
        		)
        	);
        	Sc_.set
        	(
        	        i,
        		new volScalarField
        		(
            			IOobject
            			(
               				"Sc_"+species2[i],
                			mesh.time().timeName(),
                			mesh,
                			IOobject::NO_READ,
                			IOobject::NO_WRITE
            			),
            			mesh,
            			dimless
        		)
        	);
        	SH_.set
        	(
        	        i,
        		new volScalarField
        		(
            			IOobject
            			(
               				"SH_"+species2[i],
                			mesh.time().timeName(),
                			mesh,
                			IOobject::NO_READ,
                			IOobject::AUTO_WRITE
            			),
            			mesh,
            			dimless
        		)
        	);
        }
        sh_.set
    	(
        	0,
        	new dimensionedScalar("sh_a", dimless,speciesTransportCoeffs_)
        );
        sh_.set
    	(
        	1,
        	new dimensionedScalar("sh_b", dimless,speciesTransportCoeffs_)
        );
        sh_.set
    	(
        	2,
        	new dimensionedScalar("sh_c", dimless,speciesTransportCoeffs_)
        );
}
// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::speciesTransport::correct(const int i, const volScalarField& theta_)
{

	Re_p_=max((mag(U_)+dimensionedScalar("USMALL",U_.dimensions(),SMALL))*D_pore_[0]/nuc_,ReSMALLF_);//this have to be done smarter!!
	Re_p_.correctBoundaryConditions();
		Sc_[i]=(nuc_)/D2_[i];
		Sc_[i].correctBoundaryConditions();
		SH_[i]=max(sh_[0]*pow(Re_p_,sh_[1])*pow(Sc_[i],sh_[2]),dimensionedScalar(dimless,1));
		SH_[i].correctBoundaryConditions();
		k_as_[i]=SH_[i]*(D2_[i])/(D_pore_[0]);
		k_as_[i].correctBoundaryConditions();
		C2_s_[i]=(mSTaPtr_.Psi_BV()[i]-mSTaPtr_.mDot_Wall()[i]/MW_[i])/(k_as_[i]*as_[0]*(1-theta_))+C2_[i];
		C2_s_[i].correctBoundaryConditions();
}



// ************************************************************************* //
