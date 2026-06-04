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

Class
    Foam::BearScheidegger

Description
    BearScheidegger relative velocity model

SourceFiles
    BearScheidegger.C

\*---------------------------------------------------------------------------*/

#ifndef BearScheidegger_H
#define BearScheidegger_H

#include "phaseDiffusionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseDiffusionModels
{

/*---------------------------------------------------------------------------*\
                         Class BearScheidegger Declaration
\*---------------------------------------------------------------------------*/

class BearScheidegger
:
    public phaseDiffusionModel
{
    // Private data
    
        List<word> electrodes;
        const dictionary& dict_;
        const meshObjects::gravity& g_;
        const volVectorField&   U_;
        dimensionedScalar Dp_;
        dimensionedScalar Dn_;
        //const dimensionedScalar n_;
        vector ep_; //unit vector parallel to g
        vector t_; //test vector
        vector en1_; //unit vector normal to g
        vector en2_; //unit vector normal to g
        volTensorField D_;
	    dimensionedVector DVn1_;
	    dimensionedVector DVn2_;
	    dimensionedVector DVp_;
	    
	    const scalar aL_;
	    const scalar bP_;
	    const scalar cL_;
	    const scalar dP_;
	    
	    // BSPore_ = BSPore0_*sL^nSat*U
	    const dimensionedScalar BSPore0_; //dispersivit value at 1 alphaL, in the longitudinal direction.
	    const scalar    nSat_; //BSPore0_ dependency on saturation typically negative
	    const scalar    chi_; // relation between tranverse and longitudinal longitudinal*chi=transvers dispersion, order 0.1-0.3.

        
    	volScalarField f(); //hindering fuunnnnctiiion
    	volVectorField vStokes(); //Terminal velocity model
    	
    	volTensorField UHdiff();
        volTensorField UBSdiff();

public:

    //- Runtime type information
    TypeName("BearScheidegger");


    // Constructors

        //- Construct from components
        BearScheidegger
        (
            const dictionary& dict,
            const incompressibleTwoPhaseInteractingMixture& mixture,
            const word& modelName = ""
        );


    //- Destructor
    ~BearScheidegger();


    // Member Functions

        //- Update the diffusion velocity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace phaseDiffusionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
