/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014 OpenFOAM Foundation
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

#include "mixtureViscosityModel.H"
#include "volFields.H"
#include "surfaceMesh.H"
#include "../../speciesProperties/speciesProperties.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mixtureViscosityModel, 0);
    defineRunTimeSelectionTable(mixtureViscosityModel, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mixtureViscosityModel::mixtureViscosityModel
(
    const word& name,
    const dictionary& viscosityPropertiesSub1,
    const dictionary& viscosityPropertiesSub2,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    name_(name),
    species1({"H2","O2","H2O"}),
    porousRegions({"Ne","Pe","Mem"}),
    viscosityPropertiesSub1_(viscosityPropertiesSub1),
    viscosityPropertiesSub2_(viscosityPropertiesSub2),
    muc_
    (
        "mu",dimDynamicViscosity,viscosityPropertiesSub2_
    ),
    rhoc_("rho", dimDensity, viscosityPropertiesSub2_),
    U_(U),
    phi_(phi),
    MW_(
		U_.mesh().lookupObject<speciesProperties>
        	(
            		"speciesProperties"
        	).MW()
	),
    T_(
	        U_.mesh().lookupObject<volScalarField>
        	(
            		"T"
        	)
	),
	/*mu_m
	(
    		IOobject
    		(
        		"mu_m",
        		U.mesh().time().timeName(),
        		U.mesh(),
        		IOobject::READ_IF_PRESENT,
        		IOobject::AUTO_WRITE
    		),
    		U.mesh(),
    		dimensionedScalar("mu_m__",dimDynamicViscosity,1)
	),*/
	phi1_(species1.size()*species1.size()),
	Pe_
       	(
    		U.mesh().lookupObject<volScalarField>
   	 	(
    			"Pe"
    			//(
    			//)
    		)
  	),
  	Ne_
       	(
    		U.mesh().lookupObject<volScalarField>
   	 	(
    			"Ne"
    			//(
    			//)
    		)
  	),
  	PeC_
       	(
    		U.mesh().lookupObject<volScalarField>
   	 	(
    			"PeC"
    			//(
    			//)
    		)
  	),
  	NeC_
       	(
    		U.mesh().lookupObject<volScalarField>
   	 	(
    			"NeC"
    			//(
    			//)
    		)
  	),
  	Mem_
       	(
    		U.mesh().lookupObject<volScalarField>
   	 	(
    			"Mem"
    			//(
    			//)
    		)
  	),
  	
    kr_(

	        U_.mesh().lookupObject<volScalarField>
        	(
            		"Kr"
        	)

        
        ),
        waterVapour_(viscosityPropertiesSub1_.get<bool>("waterVapour"))

	
{

}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::mixtureViscosityModel::read(const dictionary& viscosityPropertiesSub1)
{
    viscosityPropertiesSub1_ = viscosityPropertiesSub1;

    return true;
}


// ************************************************************************* //
