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

#include "capillaryD.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseDiffusionModels
{
    defineTypeNameAndDebug(capillaryD, 0);
    addToRunTimeSelectionTable(phaseDiffusionModel, capillaryD, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseDiffusionModels::capillaryD::capillaryD
(
    const dictionary& dict,
    const incompressibleTwoPhaseInteractingMixture& mixture,
    const word& modelName
)
:
    phaseDiffusionModel(dict, mixture,modelName),
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
    Solid_(
	    alphac_.mesh().lookupObject<volScalarField>
            (
            	"Solid"
            )
        )
{

}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseDiffusionModels::capillaryD::~capillaryD()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
volScalarField Foam::phaseDiffusionModels::capillaryD::Mc()
{

	return	K_*pow(alphac_,kr_)/(mixture_.muc());
}
volScalarField Foam::phaseDiffusionModels::capillaryD::Md()
{

	return	K_*pow(alphad_,kr_)/(mixture_.mud_m());
}
volScalarField Foam::phaseDiffusionModels::capillaryD::Mm()
{

	return	(mixture_.rhod()*Md()+mixture_.rhoc()*Mc())/mixture_.rho();
}
volScalarField Foam::phaseDiffusionModels::capillaryD::Ar()
{

	return	Md()*alphac_/alphad_+Mc()*alphad_/alphac_;
}
volScalarField Foam::phaseDiffusionModels::capillaryD::Cr()
{

	return	Md()/alphad_-Mc()/alphac_;
}





// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::phaseDiffusionModels::capillaryD::correct()
{
    const volScalarField& pc(mixture_.pc());
    const volScalarField& dpcds(mixture_.dpcds());
    //Info<<"min and max dpcds" << min(mixture_.dpcds())<< " : " << max(mixture_.dpcds())<<endl;
    volTensorField Drel_=
    (
        (Ar()*(-1)*dpcds-Cr()*pc)*tensor::I
    );
    volVectorField gradPcX(fvc::grad(pc) - (-1)*dpcds*fvc::grad(alphad_));
    
    BSCap_=mixture_.rhoc()/mixture_.rho()*(Md()*sqr(alphac_)+Mc()*sqr(alphad_))*gradPcX*(1-Mem_);
    Drel_.correctBoundaryConditions();
    Ddm_=(alphac_*rhoc_/mixture_.rho()) * Drel_*(1-Mem_);
    Info<<"capillaryD Ddm_ min:max  "<<endl<<"  "<<min(Ddm_)<<endl<<"  "<<max(Ddm_)<<endl;
    Info<<"capillaryD BSCap_ min:max  "<<endl<<"  "<<min(BSCap_)<<endl<<"  "<<max(BSCap_)<<endl;
    
    Info<< "max |grad(pc)|       = " << max(mag(fvc::grad(pc))) << nl;
Info<< "max |dpcdSd grad(sd)|= " << max(mag((-1)*dpcds*fvc::grad(alphad_))) << nl;
Info<< "max |gradPcX|        = " << max(mag(gradPcX)) << nl;
}


// ************************************************************************* //
