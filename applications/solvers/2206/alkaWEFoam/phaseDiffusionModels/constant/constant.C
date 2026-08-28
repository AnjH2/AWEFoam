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

#include "constant.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseDiffusionModels
{
    defineTypeNameAndDebug(constant, 0);
    addToRunTimeSelectionTable(phaseDiffusionModel, constant, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseDiffusionModels::constant::constant
(
    const dictionary& dict,
    const incompressibleTwoPhaseInteractingMixture& mixture,
    const word& modelName
)
:
    phaseDiffusionModel(dict, mixture,modelName),
    electrodes({"Ne","Pe"}),
    dict_(dict),
    g_(meshObjects::gravity::New(mixture.U().time())),
    
    rd_(electrodes.size()),
   
    D_(electrodes.size()),


    eg_("eg",(-1*g_)/mag(g_))
{
forAll(electrodes,i)
	{
	rd_.set
    	(
        	i,
        	new dimensionedScalar("rd_"+electrodes[i], dimLength,dict_)
        );
        D_.set
        (
        	i,
        	new dimensionedTensor("D_"+electrodes[i], dimless,dict_)
        );
        }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseDiffusionModels::constant::~constant()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

volVectorField Foam::phaseDiffusionModels::constant::vStokes(int j)
{
	return mag(g_)*sqr(2*rd_[j])/(18*mixture_.nucModel().nu())*eg_;
}

void Foam::phaseDiffusionModels::constant::correct()
{
    // Scale r_b |u_inf| with the prescribed anisotropy tensor and convert
    // the closure to the mixture-mass reference frame.
    Ddm_=(rhoc_/mixture_.rho())*((rd_[0]*mag(vStokes(0))*D_[0]*(Ne_*dF_+NeC_)+rd_[1]*mag(vStokes(1))*D_[1]*(Pe_*dF_+PeC_))/alphad_);

}


// ************************************************************************* //
