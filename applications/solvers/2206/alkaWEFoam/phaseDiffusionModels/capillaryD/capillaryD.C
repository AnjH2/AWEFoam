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
            	"permeabilityField"
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

	return	mixture_.rho()*K_*pow(alphac_,kr_)/(rhoc_*max(alphac_,SMALL)*mixture_.muc());
}
volScalarField Foam::phaseDiffusionModels::capillaryD::Md()
{

	return	mixture_.rho()*K_*pow(alphad_,kr_)/(rhod_*max(alphad_,SMALL)*mixture_.mud_m());
}
volScalarField Foam::phaseDiffusionModels::capillaryD::Mm()
{

	return	K_/mixture_.mu();
}







// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::phaseDiffusionModels::capillaryD::correct()
{
    const volScalarField& pc(mixture_.pc());
    const volScalarField& dpcds(mixture_.dpcds());
    Info<<"min and max dpcds" << min(mixture_.dpcds())<< " : " << max(mixture_.dpcds())<<endl;
    volTensorField Drel_=
    (
        Mm()*dpcds*tensor(1,0,0,0,1,0,0,0,1)
        +
        Mm()*pc*tensor(1,0,0,0,1,0,0,0,1)
    );
 
    Drel_.correctBoundaryConditions();
    Ddm_=(alphac_*rhoc_/mixture_.rho()) * Drel_*(1-Mem_);
    Info<<"capillaryD Ddm_ min:max  "<<endl<<"  "<<min(Ddm_)<<endl<<"  "<<max(Ddm_)<<endl;
}


// ************************************************************************* //
