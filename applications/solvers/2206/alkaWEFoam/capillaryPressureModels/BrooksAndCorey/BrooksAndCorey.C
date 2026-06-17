/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2017 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "BrooksAndCorey.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace capillaryPressureModels
{
    defineTypeNameAndDebug(BrooksAndCorey, 0);

    addToRunTimeSelectionTable
    (
        capillaryPressureModel,
        BrooksAndCorey,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::capillaryPressureModels::BrooksAndCorey::BrooksAndCorey
(
    const volScalarField& alphaWetting,
    const dictionary& dict,
    const word modelName
)
:
    capillaryPressureModel(alphaWetting,dict),
    BrooksAndCoreyCoeffsSub_(dict_.subDict(modelName + "Coeffs")),
    pc0_
    (
        IOobject
        (
            "pc0",
            alphaWetting.time().timeName(),
            alphaWetting.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),       
        alphaWetting.mesh(),
        dimensionedScalar("pc0", dimless, SMALL)
    ),
    beta_
    (
        IOobject
        (
            "beta",
            alphaWetting.time().timeName(),
            alphaWetting.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),       
        alphaWetting.mesh(),
        dimensionedScalar("beta", dimless, 1)
    ),
    alphaWetMin_("alphaWetMin",dimless,BrooksAndCoreyCoeffsSub_),
    alphaWetMax_("alphaWetMax",dimless,BrooksAndCoreyCoeffsSub_),
    alphaWetEff_
    (
        IOobject
        (
            "alphaEff",
            alphaWetting.time().timeName(),
            alphaWetting.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),       
        alphaWetting.mesh(),
        dimensionedScalar("alphaEff", dimless, 0)
    ),
    eps_
    (
    		alphaWetting.mesh().lookupObject<volScalarField>
   	 	    (
    			"eps"
    			//(
    			//)
    		)
  	),
    K_
    (
    		alphaWetting.mesh().lookupObject<volScalarField>
   	 	    (
    			"K"
    			//(
    			//)
    		)
  	),
  	Solid_
    (
    		alphaWetting.mesh().lookupObject<volScalarField>
   	 	    (
    			"Solid"
    			//(
    			//)
    		)
  	),
  	Mem_
    (
    		alphaWetting.mesh().lookupObject<volScalarField>
   	 	    (
    			"Mem"
    			//(
    			//)
    		)
  	)
{
    if (!pc0_.headerOk())
    {
    
    pc0_ =
        dimensionedScalar("pc0_Ne",dimless,BrooksAndCoreyCoeffsSub_)*Ne_
      + dimensionedScalar("pc0_Pe",dimless,BrooksAndCoreyCoeffsSub_)*Pe_
      + dimensionedScalar("pc0_Mem",dimless,BrooksAndCoreyCoeffsSub_)*Mem_
      + dimensionedScalar("pc0_min", dimless, SMALL);
        
    }
    
    if (!beta_.headerOk())
    {
    
    beta_ =
        dimensionedScalar("beta_Ne",dimless,BrooksAndCoreyCoeffsSub_)*Ne_
      + dimensionedScalar("beta_Pe",dimless,BrooksAndCoreyCoeffsSub_)*Pe_
      + dimensionedScalar("beta_Mem",dimless,BrooksAndCoreyCoeffsSub_)*Mem_
      + dimensionedScalar("beta_min", dimless, SMALL);
        
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


void Foam::capillaryPressureModels::BrooksAndCorey::correct()
{
    //! Important: We define the capillary pressure as pc = p_gas - p_water in the equations
    //! Within the equations I define pc = pgas - pliquid

   // // Calcluate effective saturation for wetting phase
   volScalarField alphaWetEffRaw_((alphaWetting_-alphaWetMin_)/(alphaWetMax_-alphaWetMin_));
   //alphaWetEff_=((alphaWetting_-alphaWetMin_)/(alphaWetMax_-alphaWetMin_));
   alphaWetEff_=min(1-SMALL,max(SMALL,alphaWetEffRaw_));
   
   //alphaWetEff_.min(1e-4);

   
   volScalarField alphaWetEffdS (
                                    (pos(alphaWetEffRaw_ - SMALL) * pos((1-SMALL) - alphaWetEffRaw_)) 
                                    / (alphaWetMax_.value()-alphaWetMin_.value())
                                );
   
   pc_=
            //(-1.0)* //ANJ
            (pc0_)*pow(alphaWetEff_, -beta_)*dimensionedScalar(dimPressure,1);
   pc_.correctBoundaryConditions();  
   dpcds_ = -1.0*
   beta_*(pc0_)*pow(alphaWetEff_, -beta_-1)*(alphaWetEffdS)*dimensionedScalar(dimPressure,1);
    dpcds_.correctBoundaryConditions(); 

}

bool Foam::capillaryPressureModels::BrooksAndCorey::read
(
    const dictionary& dict
)
{
    capillaryPressureModel::read(dict);

    BrooksAndCoreyCoeffsSub_ = dict.optionalSubDict(typeName + "Coeffs");

    //BrooksAndCoreyCoeffsSub1_.readEntry("mu", mud_);
  //  BrooksAndCoreyCoeffs_.readEntry("n", BrooksAndCoreyViscosityExponent_);
   // BrooksAndCoreyCoeffs_.readEntry("muMax", muMax_);

    return true;
}


// ************************************************************************* //
