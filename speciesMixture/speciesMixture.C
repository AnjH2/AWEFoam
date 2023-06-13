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
	//elChem_(mesh),
	/*poroM_
    	(
    		porousProperties(mesh)
    	),*/
    	species2Coeffs_(dict.optionalSubDict(phase2NamE_)),
    	species1Coeffs_(dict.optionalSubDict(phase1NamE_)),
    	species2({"H2","O2","H2O","OH"}),
    	species1({"H2","O2"}),
	z_(species2.size()),
	Ea2_(species2.size()),
	D2_ref_(species2.size()),
	T2_ref_(species2.size()),
	MW_(species2.size()+1),
	C2_(species2.size()),
	D2_(species2.size()),
	D2_eff_(species2.size()),
	C1_(species1.size()),
	//Mobility
	u_OH_("u_OH",dimensionSet ( 0, -1, 0, 0, 0, 0, 0),species2Coeffs_),
	u_K_("u_K",dimensionSet ( 0, -1, 0, 0, 0, 0, 0),species2Coeffs_),
	D1_C_star_
	(
		"C_star",
		dimless,
		species1Coeffs_
	),
	D1_a_
	(
		"Da",
		dimless,
		species1Coeffs_
	),
	D1_b_
	(
		"Db",
		dimless,
		species1Coeffs_
	),
	D1_zeta_
	(
		"zeta",
		dimless,
		species1Coeffs_
	),
	D1_ref_
	(
		"D1_ref",
		dimensionSet ( 0, 2, -1, 0, 0, 0, 0),
		species1Coeffs_
	),
    	T_
    	(
    	    IOobject
    	    (
     	       "T",
    	       mesh.time().timeName(),
     	       mesh,
     	       IOobject::MUST_READ,
    	       IOobject::NO_WRITE
    	    ),
    	    mesh 
    	),
   	 
   	epsilon_
   	(
    		mesh.lookupObject<volScalarField>
   	 	(
    			"eps"
    			//(
    			//)
    		)
  	),
    	tau_
    	(
    		mesh.lookupObject<volScalarField>
    		(
    			"tau"
    			//(
    			//)
    		)
    	),
    	alphad1_
    	(
    		mesh.lookupObject<volScalarField>
    		(
    			"alpha.gas"
    			//(
    			//)
    		)
    	)
    	/*dummy_
    	(
    	    IOobject
    	    (
     	       "___",
    	       mesh.time().timeName(),
     	       mesh,
     	       IOobject::NO_READ,
    	       IOobject::NO_WRITE
    	    ),
    	    T_/T_*0.5 
    	)
    	D1_eff_
    	(
    	    IOobject
    	    (
     	       "D1_eff",
    	       mesh.time().timeName(),
     	       mesh,
     	       IOobject::NO_READ,
    	       IOobject::AUTO_WRITE
    	    ),
    	    dummy_*D1_ref_
    	)*/
{
forAll(species2,i)
	{
		Info<< "*** Reading speciesProperties for phase2."
        	<< species2[i] << "***" << nl << endl;
       		Info<< "    Adding charges to speceis\n" << endl;
       		//z2_i="z_"+"H2";
    		z_.set
    		(
        	i,
        	new dimensionedScalar("z_"+species2[i], dimensionSet ( 0, 0, 0 , 0, 0, 0, 0),species2Coeffs_)
        	);
        	Info<< "    Adding diffusion activation constants to the speceis\n" << endl;
        	Ea2_.set
    		(
        	i,
        	new dimensionedScalar("Ea_"+species2[i], dimensionSet ( 1, 2, -2, 0, -1, 0, 0),species2Coeffs_)
        	);
        	Info<< "    Adding diffusion referennce to the speceis\n" << endl;
        	D2_ref_.set
    		(
        	i,
        	new dimensionedScalar("D_"+species2[i]+"_ref", dimensionSet ( 0, 2, -1, 0, 0, 0, 0),species2Coeffs_)
        	);
        	Info<< "    Adding temperature referennce to the speceis\n" << endl;
        	T2_ref_.set
    		(
        	i,
        	new dimensionedScalar("T_ref_"+species2[i], dimensionSet ( 0, 0, 0, 1, 0, 0, 0),species2Coeffs_)
        	);
        	Info<< "    Adding molar weight to the speceis\n" << endl;
        	MW_.set
    		(
        	i,
        	new dimensionedScalar("MW_"+species2[i], dimensionSet ( 1, 0, 0, 0, -1, 0, 0),dict)
        	);
        	Info<< "*** Reading fields for phase2"<< endl;
        	Info<< "    reading concentrations\n" << endl;
    		C2_.set
    		(
        		i,
        		new volScalarField
        		(
            			IOobject
            			(
                			"C_"+species2[i]+"."+phase2NamE_,
                			mesh.time().timeName(),
                			mesh,
                			IOobject::MUST_READ,
                			IOobject::AUTO_WRITE
            			),
            			mesh
        		)
    		);
    		Info<< "    Calculating effective diffusion\n" << endl;
    		D2_.set
    		(
        		i,
        		new volScalarField
        		(
            			IOobject
            			(
                			"D_"+species2[i]+"."+phase2NamE_,
                			mesh.time().timeName(),
                			mesh,
                			IOobject::NO_READ,
                			IOobject::NO_WRITE
            			),
            			D2_ref_[i]*exp(-Ea2_[i]/(Foam::constant::physicoChemical::R)*(1/T_-1/T2_ref_[i]))
        		)
    		);
    		D2_eff_.set
    		(
        		i,
        		new volScalarField
        		(
            			IOobject
            			(
                			"D_"+species2[i]+"_eff."+phase2NamE_,
                			mesh.time().timeName(),
                			mesh,
                			IOobject::NO_READ,
                			IOobject::NO_WRITE
            			),
            			Foam::pow(epsilon_,tau_)*Foam::pow(1-alphad1_,1-tau_)*D2_[i]
        		)
    		);
    	
	}
	Info<< "*** Reading molar weight for K***"<<" size:"<<MW_.size()<< nl << endl;
MW_.set
	(
        MW_.size()-1,
        new dimensionedScalar("MW_K", dimensionSet ( 1, 0, 0, 0, -1, 0, 0),dict)
        );
Info<< "*** Starting preparation of species 1***"<< nl << endl;

forAll(species1,i)
	{
	Info<< "*** Reading speciesProperties for phase1."
        	<< species1[i] << "***" << nl << "i="<<i<< endl;
    		C1_.set
    		(
        		i,
        		new volScalarField
        		(
            			IOobject
            			(
                			"C_"+species1[i]+"."+phase1NamE_,
                			mesh.time().timeName(),
                			mesh,
                			IOobject::MUST_READ,
                			IOobject::AUTO_WRITE
            			),
            			mesh
        		)
    		);
	}
/*Info<< "*** Calculating gasous diffusion coefficients***"<< nl << endl;
   D1_eff_=D1_ref_*(1+DeltaX1(C1_[1]/(C1_[0]+C1_[1])))/(1+DeltaX1(dummy_));
   */
    
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //




// ************************************************************************* //
