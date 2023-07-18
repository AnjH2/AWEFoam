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

#include "reactionProperties.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "fvc.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reactionProperties::reactionProperties
(
    const fvMesh& mesh
)
:
    IOdictionary
    (
        IOobject
        (
            "transportProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    	reactionCoeffs_(this->optionalSubDict("reactionCoeffs")),
    	electrodes({"Ne","Pe"}),
    	species2({"H2","O2","H2O","OH"}),
    	C_ref_(species2.size()),
	DeltaH_(2),
	H_ref_(2),
	T_H_ref_(2),
	k_H_(2),
    	psi_(species2.size()*electrodes.size()),
    	Psi_(species2.size()),
	E0_ref_(electrodes.size()),
	T0_E0_(electrodes.size()),
	Ds_(electrodes.size()),
	E0_(electrodes.size()),
	i0_ref_(electrodes.size()),
	Ea_(electrodes.size()),
	T0_i0_(electrodes.size()),
	i0_(electrodes.size()),
	alphaT_(electrodes.size()*2),
	n_
	(
		"n",
		dimensionSet( 0, 0, 0, 0, 0, 0, 0),
		reactionCoeffs_
	),
	u_K_
	(
		"u_K",
		dimensionSet( 0, -1, 0, 0, 0, 0, 0),
		reactionCoeffs_
	),
	u_OH_
	(
		"u_OH",
		dimensionSet( 0, -1, 0, 0, 0, 0, 0),
		reactionCoeffs_
	),
	tau_c_
	(
		"tau_c",
		dimensionSet( 0, 0, 1, 0, 0, 0, 0),
		reactionCoeffs_
	),
	
	
	T_
   	(
    		mesh.lookupObject<volScalarField>
   	 	(
    			"T"
    			//(
    			//)
    		)
  	),
	J_lim_
	(
		"J_lim",
		dimensionSet( 0, -2, 0, 0, 0, 1, 0),
		reactionCoeffs_
	)
{
forAll(species2,i)
	{
		Info<< "reading reference concentration for "+species2[i]+" used in butler volmer equation\n" << endl;
        	C_ref_.set
    		(
        	i,
        	new dimensionedScalar("C_"+species2[i]+"_ref", dimensionSet ( 0, -3, 0, 0, 1, 0, 0),reactionCoeffs_)
        	);
        	if (species2[i]=="H2"||species2[i]=="O2")
        	{
        		DeltaH_.set
    			(
        			i,
        			new dimensionedScalar("DeltaH_"+species2[i], dimensionSet ( 0, 0, 0, 1, 0, 0, 0),reactionCoeffs_)
        		);
        		H_ref_.set
    			(
        			i,
        			new dimensionedScalar("H_ref_"+species2[i], dimensionSet ( -1, -2, 2, 0, 1, 0, 0),reactionCoeffs_)
        		);
        		T_H_ref_.set
    			(
        			i,
        			new dimensionedScalar("T_H_ref_"+species2[i], dimensionSet ( 0, 0, 0, 1, 0, 0, 0),reactionCoeffs_)
        		);

        		k_H_.set
    			(
        			i,
        			new volScalarField
        			(
            				IOobject
            				(
               					"k_H_"+species2[i],
                				mesh.time().timeName(),
                				mesh,
                				IOobject::NO_READ,
                				IOobject::NO_WRITE
            				),
            			H_ref_[i]*exp(DeltaH_[i]*(1/T_-1/T_H_ref_[i]))
        			)
    			);
        		
        	}
        	Psi_.set
    		(
        		i,
        		new volScalarField
        		(
            			IOobject
            			(
               				"Psi_"+species2[i],
                			mesh.time().timeName(),
                			mesh,
                			IOobject::NO_READ,
                			IOobject::NO_WRITE
            			),
            			mesh,
            			dimensionSet(0,-3,-1,0,1,0,0)
        		)
    		);
	}
	
forAll(electrodes,i)
	{
		forAll(species2,j)
		{
		Info<< "reading stoichiometric coefficient for "+species2[j]+" at "+electrodes[i]+"\n" << endl;
		psi_.set
    		(
        	(i*4)+j,
        	new dimensionedScalar("psi_"+electrodes[i]+"_"+species2[j], dimless,reactionCoeffs_)
        	);
		}

	E0_ref_.set
    	(
        	i,
        	new dimensionedScalar("E0_"+electrodes[i]+"_ref", dimensionSet ( 1, 2, -3, 0, 0, -1, 0),reactionCoeffs_)
        );
        T0_E0_.set
    	(
        	i,
        	new dimensionedScalar("T0_E0_"+electrodes[i], dimensionSet ( 0, 0, 0, 1, 0, 0, 0),reactionCoeffs_)
        );
        Ds_.set
    	(
        	i,
        	new dimensionedScalar("Ds_"+electrodes[i], dimensionSet( 1, 2, -2, -1, -1, 0, 0),reactionCoeffs_)
        );
        Info<< "reading reference activation potential for "+electrodes[i]+"\n" << endl;
        E0_.set
    	(
        	i,
        	new volScalarField
        	(
            		IOobject
            		(
                		"E0_"+electrodes[i],
                		mesh.time().timeName(),
                		mesh,
                		IOobject::NO_READ,
                		IOobject::NO_WRITE
            		),
            		E0_ref_[i]+Ds_[i]/(2*Foam::constant::physicoChemical::F)*(T_-T0_E0_[i])
        	)
    	);

    	i0_ref_.set
    	(
        	i,
        	new dimensionedScalar("i0_"+electrodes[i]+"_ref", dimensionSet( 0, -2, 0, 0, 0, 1, 0),reactionCoeffs_)
        );
        T0_i0_.set
    	(
        	i,
        	new dimensionedScalar("T0_i0_"+electrodes[i], dimensionSet( 0, 0, 0, 1, 0, 0, 0),reactionCoeffs_)
        );
        Ea_.set
    	(
        	i,
        	new dimensionedScalar("Ea_"+electrodes[i], dimensionSet( 1, 2, -2, 0, -1, 0, 0),reactionCoeffs_)
        );
        i0_.set
    	(
        	i,
        	new volScalarField
        	(
            		IOobject
            		(
                		"i0_"+electrodes[i],
                		mesh.time().timeName(),
                		mesh,
                		IOobject::NO_READ,
                		IOobject::NO_WRITE
            		),
            		i0_ref_[i]*exp(-Ea_[i]/Foam::constant::physicoChemical::R*(1/T_-1/T0_i0_[i]))
        	)
    	);
    	alphaT_.set
    	(
        	i*2,
        	new dimensionedScalar("alpha_a_"+electrodes[i], dimless,reactionCoeffs_)
        );
        alphaT_.set
    	(
        	i*2+1,
        	new dimensionedScalar("alpha_c_"+electrodes[i], dimless,reactionCoeffs_)
        );
    	
	}
t_OH_=u_OH_/(u_K_+u_OH_);

}
// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //




// ************************************************************************* //
