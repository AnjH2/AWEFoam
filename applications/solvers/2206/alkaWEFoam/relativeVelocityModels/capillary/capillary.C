/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2015 OpenFOAM Foundation
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

#include "capillary.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace relativeVelocityModels
{
    defineTypeNameAndDebug(capillary, 0);
    addToRunTimeSelectionTable(relativeVelocityModel, capillary, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::relativeVelocityModels::capillary::capillary
(
    const dictionary& dict,
    const incompressibleTwoPhaseInteractingMixture& mixture,
    const word& modelName
)
:
    relativeVelocityModel(dict, mixture,modelName),
    mixture_(mixture),
    electrodes({"Ne","Pe"}),
    dict_(dict),
    
    Urel_
    (
        IOobject
        (
            modelName+"Urel",
            alphac_.time().timeName(),
            alphac_.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        alphac_.mesh(),
        dimensionedVector(dimVelocity, Zero)
    ),
    
    g_(meshObjects::gravity::New(mixture.U().time())),
   
    eps_(
	    alphac_.mesh().lookupObject<volScalarField>
            (
            	"eps"
            )
        ),
    K_(
	    alphac_.mesh().lookupObject<volScalarField>
            (
            	"permeabilityField"
            )
        ),
    kr_(mixture_.kr()),
    p_(
	    alphac_.mesh().lookupObject<volScalarField>
            (
            	"p_rgh"
            )
        ),
    Solid_(
	    alphac_.mesh().lookupObject<volScalarField>
            (
            	"Solid"
            )
        )
{

}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::relativeVelocityModels::capillary::~capillary()
{}

// * * * * * * * * * * * * * * Private Functions  * * * * * * * * * * * * * * //
volScalarField Foam::relativeVelocityModels::capillary::Mc()
{

	return	mixture_.rho()*K_*pow(alphac_,kr_)/(rhoc_*max(alphac_,SMALL)*mixture_.muc());
}
volScalarField Foam::relativeVelocityModels::capillary::Md()
{

	return	mixture_.rho()*K_*pow(alphad_,kr_)/(rhod_*max(alphad_,SMALL)*mixture_.mud_m());
}
volScalarField Foam::relativeVelocityModels::capillary::Mm()
{

	return	K_/mixture_.mu();
}



// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::relativeVelocityModels::capillary::correctRelativeVelocity()
{

    const volScalarField& pc(mixture_.pc());
    
    //volVectorField& Urel = this->Urel_;

    Urel_=
    (//(-Md()+Mc())*fvc::grad(p_)+//it is simply to strong!
    (-rhoc_*Mc()+rhod_*Md())*g_
    //-(Md()*alphac_+Mc()*alphad_)*fvc::grad(pc)
    //+(Mc()-Md())*pc*fvc::grad(alphac_)
    )*Solid_*-1;
 
    Urel_.correctBoundaryConditions();
}

void Foam::relativeVelocityModels::capillary::correct()
{    
    dModel_->correct();
    correctRelativeVelocity();
    Udm_ = (alphac_*rhoc_/mixture_.rho()) * Urel_*(1-hF_)*(1-Mem_+VSMALL);
    

    /*
    Udm_.component(0) = pow(2,0.5)*pow((sigma_*g_.component(0)*(rhod_-rhoc_))/(pow(rhoc_,2)),0.25)*pow(1-alphad_,1.75);
    Info<<Udm_.component(0)<<endl;
    Udm_.component(1) = pow(2,0.5)*pow((sigma_*g_.component(1)*(rhod_-rhoc_))/(pow(rhoc_,2)),0.25)*pow(1-alphad_,1.75);
    Info<<Udm_.component(1)<<endl;
    Udm_.component(2) = pow(2,0.5)*pow((sigma_*g_.component(2)*(rhod_-rhoc_))/(pow(rhoc_,2)),0.25)*pow(1-alphad_,1.75);
    Info<<Udm_.component(2)<<endl;
    Info<<((sigma_*g_*(rhoc_-rhod_))/(pow(rhoc_,2))).component(0)<<endl;*/
    
    //Ddm_=(rhoc_/rho())*(rd_[0]*mag(vcapillary(0))*D_[0]*(Ne_*dF_+NeC_)+rd_[1]*mag(vcapillary(1))*D_[1]*(Pe_*dF_+PeC_))*f()/alphad_;
    //Ddm_.correctBoundaryConditions();
}


// ************************************************************************* //
