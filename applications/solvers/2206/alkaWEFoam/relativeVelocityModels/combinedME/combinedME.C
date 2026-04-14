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

#include "combinedME.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace relativeVelocityModels
{
    defineTypeNameAndDebug(combinedME, 0);
    addToRunTimeSelectionTable(relativeVelocityModel, combinedME, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::relativeVelocityModels::combinedME::combinedME
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
    g_(meshObjects::gravity::New(mixture.U().time())),
    
    rd_(electrodes.size()),
	//(
	//	"R_DB",
	//	dimensionSet ( 0, 1, 0, 0, 0, 0,0),
	//	dict
	//),
    n_("n",dimless,dict_),
    rhoc_(mixture.rhoc()),
    rhod_(mixture.rhod()),
      
    eps_(
	    alphac_.mesh().lookupObject<volScalarField>
            (
            	"eps"
            )
        ),

    eg_("eg",(-1*g_)/mag(g_)),
    C_("C",dimless,dict_),
    alphaStar_("alphaStar",dimless,dict_)
{
forAll(electrodes,i)
	{
	rd_.set
    	(
        	i,
        	new dimensionedScalar("rd_"+electrodes[i], dimLength,dict_)
        );
        }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::relativeVelocityModels::combinedME::~combinedME()
{}

// * * * * * * * * * * * * * * Private Functions  * * * * * * * * * * * * * * //
volScalarField Foam::relativeVelocityModels::combinedME::f()
{
	
	return pow(1-alphad_,n_);
}

dimensionedVector Foam::relativeVelocityModels::combinedME::vStokes(int j)
{
	
	return mag(g_)*sqr(2*rd_[j])/(18*mixture_.nuc())*eg_;
}



volVectorField Foam::relativeVelocityModels::combinedME::UStokes(int j)
{
	
	return f()*mag(vStokes(j))*(1-(Ne_+Pe_)*hF_)*eg_;
}
volTensorField Foam::relativeVelocityModels::combinedME::Usp(int j)
{
	
	return (pow(2*rd_[j],2)/(18*alphad_*mixture_.muc()))*exp(-C_*(alphaStar_-alphad_))*dimensionedTensor(dimPressure, One);
}






// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::relativeVelocityModels::combinedME::correct()
{    
    dModel_->correct();
    Udm_ = (1-Mem_+VSMALL)*(rhoc_/rho())*((UStokes(0))*(Ne_+NeC_)+(UStokes(1))*(Pe_+PeC_));
    

    /*
    Udm_.component(0) = pow(2,0.5)*pow((sigma_*g_.component(0)*(rhod_-rhoc_))/(pow(rhoc_,2)),0.25)*pow(1-alphad_,1.75);
    Info<<Udm_.component(0)<<endl;
    Udm_.component(1) = pow(2,0.5)*pow((sigma_*g_.component(1)*(rhod_-rhoc_))/(pow(rhoc_,2)),0.25)*pow(1-alphad_,1.75);
    Info<<Udm_.component(1)<<endl;
    Udm_.component(2) = pow(2,0.5)*pow((sigma_*g_.component(2)*(rhod_-rhoc_))/(pow(rhoc_,2)),0.25)*pow(1-alphad_,1.75);
    Info<<Udm_.component(2)<<endl;
    Info<<((sigma_*g_*(rhoc_-rhod_))/(pow(rhoc_,2))).component(0)<<endl;*/
    
    //Ddm_=(rhoc_/rho())*(((rd_[0]*mag(vStokes(0))*D_[0]+Usp(0))*(Ne_*dF_+NeC_)+(rd_[1]*mag(vStokes(1))*D_[1])*(Pe_*dF_+PeC_))*f()/alphad_+(Usp(1)*(Pe_*dF_+PeC_)+Usp(0)*(Ne_*dF_+NeC_)));
    //Ddm_.correctBoundaryConditions();
}


// ************************************************************************* //
