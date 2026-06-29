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

#include "capillaryDiscontinuous.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace relativeVelocityModels
{
    defineTypeNameAndDebug(capillaryDiscontinuous, 0);
    addToRunTimeSelectionTable(relativeVelocityModel, capillaryDiscontinuous, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::relativeVelocityModels::capillaryDiscontinuous::capillaryDiscontinuous
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
    CBeta_("CBeta",dimless,dict_),
    
    g_(meshObjects::gravity::New(mixture_.U().time())),
   
    eps_(
	    alphac_.mesh().lookupObject<volScalarField>
            (
            	"eps"
            )
        ),
    K_(
	    alphac_.mesh().lookupObject<volScalarField>
            (
            	"K"
            )
        ),
    kr_(mixture_.kr()),
    p_(
	    alphac_.mesh().lookupObject<volScalarField>
            (
            	"p"
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

Foam::relativeVelocityModels::capillaryDiscontinuous::~capillaryDiscontinuous()
{}

// * * * * * * * * * * * * * * Private Functions  * * * * * * * * * * * * * * //
volScalarField Foam::relativeVelocityModels::capillaryDiscontinuous::Mc()
{

	return	K_*pow(alphac_,kr_)/(mixture_.muc());
}
volScalarField Foam::relativeVelocityModels::capillaryDiscontinuous::Md()
{

	return	K_*pow(alphad_,kr_)/(mixture_.mud_m());
}
volScalarField Foam::relativeVelocityModels::capillaryDiscontinuous::Mm()
{

	return	(mixture_.rhoc()*Mc())/mixture_.rho();
}

volScalarField Foam::relativeVelocityModels::capillaryDiscontinuous::rhoS()
{
    return (pow(mixture_.rhod(),2)*Md()+pow(mixture_.rhoc(),2)*Mc())
            / (mixture_.rhod()*Md()+mixture_.rhoc()*Mc());
}
volScalarField Foam::relativeVelocityModels::capillaryDiscontinuous::wd()
{
    return mixture_.rhod()*Md()/ (mixture_.rhod()*Md()+mixture_.rhoc()*Mc());
}
volScalarField Foam::relativeVelocityModels::capillaryDiscontinuous::wc()
{
    return mixture_.rhoc()*Mc()/ (mixture_.rhod()*Md()+mixture_.rhoc()*Mc());
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::relativeVelocityModels::capillaryDiscontinuous::correctRelativeVelocity()
{
Info<<"re1"<<endl;
        const volScalarField& pc = mixture_.pc();
        volScalarField beta0=CBeta_/sqrt(K_)*Solid_;
        volScalarField betaC=beta0/pow(max(alphac_, 1e-4), 3);
        volScalarField betaD=beta0/pow(max(alphad_, 1e-4), 3);


Info<<"re2"<<endl;
        // Initial Darcy phase fluxes
        volVectorField qd
        (
            "qd",
            Md()*Solid_
           *(
              - fvc::grad(alphac_*pc,"grad(pc)")
            )
        );


Info<<"re3"<<endl;
Info<<"max(mag(qd)) before Forch:   "<<max(mag(qd))<<endl;
        // Forchheimer magnitude correction
        qd *=
            2.0
           /(
                1.0
              + sqrt
                (
                    1.0
                  + 4.0*mixture_.rhod()*betaD*Md()*mag(qd)
                )
            );


Info<<"re4"<<endl;
Info<<"max(mag(qd)) after Forch:   "<<max(mag(qd))<<endl;
        // Porosity-scaled relative velocity
        Urel_ =
              qd/max(alphad_, 1e-4);
        volVectorField Gd
        (
            "Gd",
            fvc::grad(alphac_*pc,"grad(pc)")
        );


Info<<"re5"<<endl;
    
    Info<< "max |grad(pd)|     = "
    << max(mag(fvc::grad(alphac_*pc,"grad(pc)"))*Solid_) << nl;

    Info<< "max |rho_d*g|      = "
        << max(mag(mixture_.rhod()*g_)*Solid_) << nl;

    Info<< "max |Gd|           = "
        << max(mag(Gd)*Solid_) << nl;

    Info<< "max |grad(pcPhase)|= "
        << max(mag(fvc::grad( - alphad_*pc,"grad(pc)"))*Solid_) << nl;



        

Info<<"re6"<<endl;
volScalarField ReKd
(
    "ReKd",
    mixture_.rhod()*mag(qd)*sqrt(K_*(1-Mem_)*Solid_)/mixture_.mud_m()
);



Info<< "ReKd mean/max = "
    << gAverage(ReKd) << "  " << max(ReKd*Solid_) << nl;

  
    
const vector gHat = g_.value()/mag(g_.value());




    
}

void Foam::relativeVelocityModels::capillaryDiscontinuous::correct()
{    
    
    dModel_->correct();
    
    correctRelativeVelocity();
    
    Udm_ = (alphac_*rhoc_/mixture_.rho()) * Urel_;
    
    //F_=((rhoS()-mixture_.rho())*g_-wd()*fvc::grad(alphac_*mixture_.pc(),"grad(pc)")+wc()*fvc::grad(alphad_*mixture_.pc(),"grad(pc)"));
    F_=(
            // (rhoStar - rho)*g . Sf
            //  fvc::interpolate(rhoS() - mixture_.rho())
            // *(g_ & alphad_.mesh().Sf())

            // -wd*grad(alphac*pc) . Sf
            - fvc::interpolate(wd())
             *fvc::snGrad(alphac_*mixture_.pc())
             *alphad_.mesh().magSf()

            // +wc*grad(alphad*pc) . Sf
          //  + fvc::interpolate(wc())
          //   *fvc::snGrad(alphad_*mixture_.pc())
          //   *alphad_.mesh().magSf()
    );

}


// ************************************************************************* //
