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

#include "speciesProperties.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "fvc.H"




// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::speciesProperties::speciesProperties
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    IOdictionary
    (
        IOobject
        	(
            	"speciesProperties",
            	mesh.time().constant(),
            	mesh,
            	IOobject::MUST_READ_IF_MODIFIED,
            	IOobject::NO_WRITE
        	)
    	),
	phase1NamE_(dict.get<wordList>("phases")[0]),
    	phase2NamE_(dict.get<wordList>("phases")[1]),
    	species2Coeffs_(this->optionalSubDict(phase2NamE_)),
    	species1Coeffs_(this->optionalSubDict(phase1NamE_)),
    	species2({"H2","O2","H2O","OH","K"}),
    	species1({"H2","O2","H2O"}),
	z_(species2.size()),
	Ea2_(species2.size()),
	D2_ref_(species2.size()),
	T2_ref_(species2.size()),
	MW_(species2.size()),
	C2_(species2.size()),
	D2_(species2.size()),
	D2_eff_(species2.size()),
	t_(species2.size()),
	D2_ambi_
    (
    	    IOobject
    	    (
     	       "D2_ambi",
    	       mesh.time().timeName(),
     	       mesh,
     	       IOobject::NO_READ,
    	       IOobject::AUTO_WRITE
    	    ),
    	    mesh,
    	    dimensionedScalar(dimVelocity*dimLength,Zero)
    ),
    C2_T_
    (
    	    IOobject
    	    (
     	       "C2_T",
    	       mesh.time().timeName(),
     	       mesh,
     	       IOobject::NO_READ,
    	       IOobject::NO_WRITE
    	    ),
    	    mesh,
    	    dimensionedScalar(dimMoles/dimVolume,Zero)
    ),
    C2_0_
    (
    	    IOobject
    	    (
     	       "C2_0",
    	       mesh.time().timeName(),
     	       mesh,
     	       IOobject::NO_READ,
    	       IOobject::AUTO_WRITE
    	    ),
    	    mesh,
    	    dimensionedScalar(dimMoles/dimVolume,Zero)
    ),
	//Mobility
	u_OH_("u_OH",dimensionSet ( 0, -1, 0, 0, 0, 0, 0),species2Coeffs_),
	u_K_("u_K",dimensionSet ( 0, -1, 0, 0, 0, 0, 0),species2Coeffs_),

    	T_
    	(
    	    IOobject
    	    (
     	       "T",
    	       mesh.time().timeName(),
     	       mesh,
     	       IOobject::MUST_READ,
    	       IOobject::AUTO_WRITE
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
    	Mem_
    	(
    		mesh.lookupObject<volScalarField>
    		(
    			"Mem"
    			//(
    			//)
    		)
    	)

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
            			Foam::pow(epsilon_,tau_)*D2_[i]
        		)
    		);
    }
    forAll(species2,i)
	{
    	if (i==3){
    	    t_.set
    		(
        		i,
        		new volScalarField
        		(
            			IOobject
            			(
                			"t_"+species2[i],
                			mesh.time().timeName(),
                			mesh,
                			IOobject::NO_READ,
                			IOobject::NO_WRITE
            			),
            			-1*z_[i]*D2_[i]/(z_[4]*D2_[4]-z_[3]*D2_[3])
        		)
    		);
    	
    	}else if (i==4){
    	    t_.set
    		(
        		i,
        		new volScalarField
        		(
            			IOobject
            			(
                			"t_"+species2[i],
                			mesh.time().timeName(),
                			mesh,
                			IOobject::NO_READ,
                			IOobject::NO_WRITE
            			),
            			mag(z_[i])*D2_[i]/(z_[4]*D2_[4]-z_[3]*D2_[3])
        		)
    		);
    	
    	}
	}
/*	Info<< "*** Reading molar weight for K***"<<" size:"<<MW_.size()<< nl << endl;
MW_.set
	(
        MW_.size()-1,
        new dimensionedScalar("MW_K", dimensionSet ( 1, 0, 0, 0, -1, 0, 0),dict)
        );
*/
D2_ambi_=(z_[4]-z_[3])*D2_[3]*D2_[4]/(-1*z_[3]*D2_[3]+z_[4]*D2_[4]);

}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
volScalarField Foam::speciesProperties::tp_Merenkov()
{
    return 0.26-0.047*(sqrt(C2_[3]/dimensionedScalar(dimMoles/dimVolume,1000))+1);
}


void Foam::speciesProperties::correct()
{
    C2_[4]=C2_[3];
    C2_T_=dimensionedScalar(dimMoles/dimVolume,Zero);
    C2_0_=dimensionedScalar(dimMoles/dimVolume,Zero);
    forAll(species2,i)
    {
        C2_T_+=C2_[i];
        if (i<3)
        {
            C2_0_+=C2_[i];
        }
    }    
}



// ************************************************************************* //
