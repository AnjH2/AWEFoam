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

Foam::relativeVelocityModels::capillary::~capillary()
{}

// * * * * * * * * * * * * * * Private Functions  * * * * * * * * * * * * * * //
volScalarField Foam::relativeVelocityModels::capillary::Mc()
{

	return	K_*pow(alphac_,kr_)/(mixture_.muc());
}
volScalarField Foam::relativeVelocityModels::capillary::Md()
{

	return	K_*pow(alphad_,kr_)/(mixture_.mud_m());
}
volScalarField Foam::relativeVelocityModels::capillary::Mm()
{

	return	(mixture_.rhod()*Md()+mixture_.rhoc()*Mc())/mixture_.rho();
}

volScalarField Foam::relativeVelocityModels::capillary::rhoS()
{
    return (pow(mixture_.rhod(),2)*Md()+pow(mixture_.rhoc(),2)*Mc())
            / (mixture_.rhod()*Md()+mixture_.rhoc()*Mc());
}
volScalarField Foam::relativeVelocityModels::capillary::wd()
{
    return mixture_.rhod()*Md()/ (mixture_.rhod()*Md()+mixture_.rhoc()*Mc());
}
volScalarField Foam::relativeVelocityModels::capillary::wc()
{
    return mixture_.rhoc()*Mc()/ (mixture_.rhod()*Md()+mixture_.rhoc()*Mc());
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::relativeVelocityModels::capillary::correctRelativeVelocity()
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
                mixture_.rhod()*g_
              - fvc::grad( alphac_*pc,"grad(pc)")
            )
        );

        volVectorField qc
        (
            "qc",
            Mc()*Solid_
           *(
                mixture_.rhoc()*g_
              - fvc::grad(- alphad_*pc,"grad(pc)")
            )
        );
Info<<"re3"<<endl;
Info<<"max(mag(qd)) before Forch:   "<<max(mag(qd))<<endl;
Info<<"max(mag(qc)) before Forch:   "<<max(mag(qc))<<endl;
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

        qc *=
            2.0
           /(
                1.0
              + sqrt
                (
                    1.0
                  + 4.0*mixture_.rhoc()*betaC*Mc()*mag(qc)
                )
            );
Info<<"re4"<<endl;
Info<<"max(mag(qd)) after Forch:   "<<max(mag(qd))<<endl;
Info<<"max(mag(qc)) after Forch:   "<<max(mag(qc))<<endl;
        // Porosity-scaled relative velocity
        Urel_ =
              qd/max(alphad_, 1e-4)
            - qc/max(alphac_, 1e-4);
        volVectorField Gd
        (
            "Gd",
            fvc::grad( alphac_*pc,"grad(pc)")
          - mixture_.rhod()*g_
        );

        volVectorField Gc
        (
            "Gc",
            fvc::grad( - alphad_*pc,"grad(pc)")
          - mixture_.rhoc()*g_
        );
Info<<"re5"<<endl;
    
    Info<< "max |grad(pd)|     = "
    << max(mag(fvc::grad( alphac_*pc,"grad(pc)"))*(Solid_*(1-Mem_))) << nl;

    Info<< "max |rho_d*g|      = "
        << max(mag(mixture_.rhod()*g_)*(Solid_*(1-Mem_))) << nl;

    Info<< "max |Gd|           = "
        << max(mag(Gd)*(Solid_*(1-Mem_))) << nl;

    Info<< "max |grad(pcPhase)|= "
        << max(mag(fvc::grad( - alphad_*pc,"grad(pc)"))*(Solid_*(1-Mem_))) << nl;

    Info<< "max |rho_c*g|      = "
        << mag(mixture_.rhoc()*g_) << nl;

    Info<< "max |Gc|           = "
        << max(mag(Gc)*(Solid_*(1-Mem_))) << nl;
        

Info<<"re6"<<endl;
volScalarField ReKd
(
    "ReKd",
    mixture_.rhod()*mag(qd)*sqrt(K_*(1-Mem_)*(Solid_*(1-Mem_)))/mixture_.mud_m()
);

volScalarField ReKc
(
    "ReKc",
    mixture_.rhoc()*mag(qc)*sqrt(K_*(1-Mem_)*(Solid_*(1-Mem_)))/mixture_.muc()
);

Info<< "ReKd mean/max = "
    << gAverage(ReKd) << "  " << max(ReKd*(Solid_*(1-Mem_))) << nl;

Info<< "ReKc mean/max = "
    << gAverage(ReKc) << "  " << max(ReKc*(Solid_*(1-Mem_))) << nl;
    
    
const vector gHat = g_.value()/mag(g_.value());

volScalarField GcAlongG
(
    "GcAlongG",
    ((mixture_.rhoc()*g_ - fvc::grad( - alphad_*pc)) & gHat)
   *(Solid_*(1-Mem_))
);

Info<< "Gc along g min/max = "
    << min(GcAlongG.primitiveField()) << "  "
    << max(GcAlongG.primitiveField()) << nl;



volVectorField GcPerpendicular
(
    "GcPerpendicular",
    Gc - (Gc & gHat)*gHat
);

Info<< "max |Gc along g| = "
    << max(mag(Gc & gHat)().primitiveField()) << nl;

Info<< "max |Gc perpendicular| = "
    << max((mag(GcPerpendicular)*(Solid_*(1-Mem_)))().primitiveField()) << nl; 
    
    
}

void Foam::relativeVelocityModels::capillary::correct()
{    
    
    dModel_->correct();
    
    correctRelativeVelocity();
    
    Udm_ = (alphac_*rhoc_/mixture_.rho()) * Urel_;
    
    //F_=((rhoS()-mixture_.rho())*g_-wd()*fvc::grad(alphac_*mixture_.pc(),"grad(pc)")+wc()*fvc::grad(alphad_*mixture_.pc(),"grad(pc)"));
    F_=(
            // (rhoStar - rho)*g . Sf
              fvc::interpolate(rhoS() - mixture_.rho())
             *(g_ & alphad_.mesh().Sf())

            // -wd*grad(alphac*pc) . Sf
            - fvc::interpolate(wd())
             *fvc::snGrad(alphac_*mixture_.pc())
             *alphad_.mesh().magSf()

            // +wc*grad(alphad*pc) . Sf
            + fvc::interpolate(wc())
             *fvc::snGrad(alphad_*mixture_.pc())
             *alphad_.mesh().magSf()
    );

}


// ************************************************************************* //
