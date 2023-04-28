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

#include "speciesMixture.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::speciesMixture::speciesMixture
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phase1NamE_,
    const word& phase2NamE_
)
:
    	speciesCoeffs_(dict.optionalSubDict(phase1NamE_)),
    	
    	z_OH_("z_OH",dimensionSet ( 0, 0, 0 , 0, 0, 0, 0),speciesCoeffs_),
    	z_H2O_("z_H2O",dimensionSet ( 0, 0, 0 , 0, 0, 0, 0),speciesCoeffs_),
    	z_H2_("z_H2",dimensionSet ( 0, 0, 0 , 0, 0, 0, 0),speciesCoeffs_),
    	z_O2_("z_O2",dimensionSet ( 0, 0, 0 , 0, 0, 0, 0),speciesCoeffs_),
	//Mobility
	u_OH_("u_OH",dimensionSet ( 0, -1, 0, 0, 0, 0, 0),speciesCoeffs_),
	u_K_("u_K",dimensionSet ( 0, -1, 0, 0, 0, 0, 0),speciesCoeffs_),
	
	//diffusion related constants
	Ea_H2O_("Ea_H2O",dimensionSet ( 1, 2, -2, 0, -1, 0, 0),speciesCoeffs_),
	Ea_OH_("Ea_OH",dimensionSet ( 1, 2, -2, 0, -1, 0, 0),speciesCoeffs_),
	Ea_H2_("Ea_H2",dimensionSet ( 1, 2, -2, 0, -1, 0, 0),speciesCoeffs_),
	Ea_O2_("Ea_O2",dimensionSet ( 1, 2, -2, 0, -1, 0, 0),speciesCoeffs_),

	D_H2O_ref_("D_H2O_ref",dimensionSet ( 0, 2, -1, 0, 0, 0, 0),speciesCoeffs_),
	D_OH_ref_("D_OH_ref",dimensionSet ( 0, 2, -1, 0, 0, 0, 0),speciesCoeffs_),
	D_H2_ref_("D_H2_ref",dimensionSet ( 0, 2, -1, 0, 0, 0, 0),speciesCoeffs_),
	D_O2_ref_("D_O2_ref",dimensionSet ( 0, 2, -1, 0, 0, 0, 0),speciesCoeffs_),

	T_ref_H2_("T_ref_H2",dimensionSet ( 0, 0, 0, 1, 0, 0, 0),speciesCoeffs_),
	T_ref_O2_("T_ref_O2",dimensionSet ( 0, 0, 0, 1, 0, 0, 0),speciesCoeffs_),
	T_ref_H2O_("T_ref_H2O",dimensionSet ( 0, 0, 0, 1, 0, 0, 0),speciesCoeffs_),
	T_ref_OH_("T_ref_OH",dimensionSet ( 0, 0, 0, 1, 0, 0, 0),speciesCoeffs_),
	
	//mole weight
	MW_H2_("MW_H2",dimensionSet ( 1, 0, 0, 0, -1, 0, 0),speciesCoeffs_),
	MW_O2_("MW_O2",dimensionSet ( 1, 0, 0, 0, -1, 0, 0),speciesCoeffs_),
	MW_OH_("MW_OH",dimensionSet ( 1, 0, 0, 0, -1, 0, 0),speciesCoeffs_),
	MW_H2O_("MW_H2O",dimensionSet ( 1, 0, 0, 0, -1, 0, 0),speciesCoeffs_),
	MW_K_("MW_K",dimensionSet ( 1, 0, 0, 0, -1, 0, 0),speciesCoeffs_),


    C_H2_1_
    (
        IOobject
        (
            IOobject::groupName("C_H2_", phase1NamE_),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    
    C_O2_1_
    (
        IOobject
        (
            IOobject::groupName("C_O2_", phase1NamE_),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    
    C_OH_1_
    (
        IOobject
        (
            IOobject::groupName("C_OH_", phase1NamE_),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    
    C_H2O_1_
    (
        IOobject
        (
            IOobject::groupName("C_H2O_", phase1NamE_),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    C_H2_2_
    (
        IOobject
        (
            IOobject::groupName("C_H2_", phase2NamE_),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    
    C_O2_2_
    (
        IOobject
        (
            IOobject::groupName("C_O2_", phase2NamE_),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
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
    
    epsilon_
    (
    	mesh.lookupObject<volScalarField>
    	(
    		"epsilon_"
    		//(
    		//)
    	)
    ),
    tau_
    (
    	mesh.lookupObject<volScalarField>
    	(
    		"tau_"
    		//(
    		//)
    	)
    ),
    //const volScalarField& T = mesh.lookupObject<volScalarField>("T"),
    //const volScalarField& epsilon_ = mesh.lookupObject<volScalarField>("epsilon"),
    //const volScalarField& tau_ = mesh.lookupObject<volScalarField>("tau"),
    //const dimensionedScalar& R_=Foam::constant::physicoChemical::R,
    //const dimensionedScalar& F_ = Foam::constant::physicoChemical::F;
    
	D_H2O
	(
	    IOobject
 	   (
	        "D_H2O",
	        mesh.time().timeName(),
	        mesh,
	        IOobject::NO_READ,
	        IOobject::NO_WRITE
	    ),
	    D_H2O_ref_*exp(-Ea_H2O_/(Foam::constant::physicoChemical::R)*(1/T_-1/T_ref_H2O_))
	),
	D_OH
	(
	    IOobject
 	   (
	        "D_OH",
	        mesh.time().timeName(),
	        mesh,
	        IOobject::NO_READ,
	        IOobject::NO_WRITE
	    ),
	    D_OH_ref_*exp(-Ea_OH_/(Foam::constant::physicoChemical::R)*(1/T_-1/T_ref_OH_))
	),
	D_H2
	(
	    IOobject
 	   (
	        "D_H2",
	        mesh.time().timeName(),
	        mesh,
	        IOobject::NO_READ,
	        IOobject::NO_WRITE
	    ),
	    D_H2_ref_*exp(-Ea_H2_/(Foam::constant::physicoChemical::R)*(1/T_-1/T_ref_H2_))
	),
	D_O2
	(
	    IOobject
 	   (
	        "D_O2",
	        mesh.time().timeName(),
	        mesh,
	        IOobject::NO_READ,
	        IOobject::NO_WRITE
	    ),
	    D_O2_ref_*exp(-Ea_O2_/(Foam::constant::physicoChemical::R)*(1/T_-1/T_ref_O2_))
	),
	//calculating effective difusivity
	D_H2O_eff
	(
	    IOobject
 	   (
	        "D_H2O_eff",
	        mesh.time().timeName(),
	        mesh,
	        IOobject::NO_READ,
	        IOobject::NO_WRITE
	    ),
	    Foam::pow(epsilon_,tau_)*Foam::pow(1,1-tau_)*D_H2O
	),
	D_OH_eff
	(
	    IOobject
 	   (
	        "D_OH_eff",
	        mesh.time().timeName(),
	        mesh,
	        IOobject::NO_READ,
	        IOobject::NO_WRITE
	    ),
	    Foam::pow(epsilon_,tau_)*Foam::pow(1,1-tau_)*D_OH
	),
	D_H2_eff
	(
	    IOobject
 	   (
	        "D_H2_eff",
	        mesh.time().timeName(),
	        mesh,
	        IOobject::NO_READ,
	        IOobject::NO_WRITE
	    ),
	    Foam::pow(epsilon_,tau_)*Foam::pow(1,1-tau_)*D_H2
	),
	D_O2_eff
	(
	    IOobject
 	   (
	        "D_O2_eff",
	        mesh.time().timeName(),
	        mesh,
	        IOobject::NO_READ,
	        IOobject::NO_WRITE
	    ),
	    Foam::pow(epsilon_,tau_)*Foam::pow(1,1-tau_)*D_O2
	)
    
    
    
{}
// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //




// ************************************************************************* //
