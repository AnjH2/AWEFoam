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

#include "leverett.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace capillaryPressureModels
{
    defineTypeNameAndDebug(leverett, 0);

    addToRunTimeSelectionTable
    (
        capillaryPressureModel,
        leverett,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::capillaryPressureModels::leverett::leverett
(
    const volScalarField& alphaWetting,
    const dictionary& dict,
    const word modelName
)
:
    capillaryPressureModel(alphaWetting,dict),
    leverettCoeffsSub_(dict_.subDict(modelName + "Coeffs")),
    sigma_("sigma",dimEnergy/dimArea,leverettCoeffsSub_),
    theta_(leverettCoeffsSub_.get<scalar>("thetaDeg")),
    alphaWetMin_("alphaWetMin",dimless,leverettCoeffsSub_),
    alphaWetMax_("alphaWetMax",dimless,leverettCoeffsSub_),
    alphaWetEff_
    (
        IOobject
        (
            "alphaEff",
            alphaWetting.time().timeName(),
            alphaWetting.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
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
    			"permeabilityField"
    			//(
    			//)
    		)
  	)
{

}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


void Foam::capillaryPressureModels::leverett::correct()
{
    //! Important: We define the capillary pressure as pc = p_gas - p_water in the equations
    //! Within the equations I define pc = pgas - pliquid

   // // Calcluate effective saturation for wetting phase
   alphaWetEff_=(alphaWetting_-alphaWetMin_)/(alphaWetMax_-alphaWetMin_);
   alphaWetEff_.max(1e-4);
   alphaWetEff_.min(1-1e-4);
   scalar thetaRad = theta_ * (Foam::constant::mathematical::pi/ 180);
    
   volScalarField sqrtK(sqrt(1/K_));
   pc_=
            //(-1)* ANJ
            sigma_*cos(thetaRad)*sqrt(eps_)*sqrtK
           *(
              thetaRad < Foam::constant::mathematical::pi/2.0
               ? (1.42*(1.0 - alphaWetEff_) - 2.12*pow(1.0 - alphaWetEff_, 2) + 1.26*pow(1.0 - alphaWetEff_, 3))
               : (1.42*alphaWetEff_ - 2.12*pow(alphaWetEff_, 2) + 1.26*pow(alphaWetEff_, 3))
              // //! Change to pc = pl - pg
              // ? (1.42*alphaWater[cellI] - 2.12*pow(alphaWater[cellI], 2) + 1.26*pow(alphaWater[cellI], 3))
              // : (1.42*(1.0 - alphaWater[cellI]) - 2.12*pow(1.0 - alphaWater[cellI], 2) + 1.26*pow(1.0 - alphaWater[cellI], 3))
           );
   pc_.correctBoundaryConditions();  
   dpcds_ =
           //! Change to pc = pl - pg  (For IMPES)
          //(-1) * ANJ
           sigma_*cos(thetaRad)*sqrt(eps_)*sqrtK
           *(
               thetaRad < Foam::constant::mathematical::pi/2.0
               ? (-1.42 + 2.0*2.12*(1.0 - alphaWetEff_) - 3.0*1.26*pow(1.0 - alphaWetEff_, 2))
               : (1.42 -2.0*2.12*alphaWetEff_ + 3.0*1.26*pow(alphaWetEff_, 2))
              // ? (1.42 -2.0*2.12*alphaWater[cellI] + 3.0*1.26*pow(alphaWater[cellI], 2))
              // : (-1.42 + 2.0*2.12*(1.0 - alphaWater[cellI]) - 3.0*1.26*pow(1.0 - alphaWater[cellI], 2))
           );
    dpcds_.correctBoundaryConditions(); 

}

bool Foam::capillaryPressureModels::leverett::read
(
    const dictionary& dict
)
{
    capillaryPressureModel::read(dict);

    leverettCoeffsSub_ = dict.optionalSubDict(typeName + "Coeffs");

    //leverettCoeffsSub1_.readEntry("mu", mud_);
  //  leverettCoeffs_.readEntry("n", leverettViscosityExponent_);
   // leverettCoeffs_.readEntry("muMax", muMax_);

    return true;
}


// ************************************************************************* //
