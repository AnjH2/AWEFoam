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
    Foam::simpleP

Description
    simpleP relative velocity model

SourceFiles
    simpleP.C

\*---------------------------------------------------------------------------*/

#ifndef simpleP_H
#define simpleP_H

#include "phaseDiffusionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseDiffusionModels
{

/*---------------------------------------------------------------------------*\
                         Class simpleP Declaration
\*---------------------------------------------------------------------------*/

class simpleP
:
    public phaseDiffusionModel
{
    // Private data
    
        
        const meshObjects::gravity& g_;
        dimensionedScalar Dp_;
        dimensionedScalar Dn_;
        //const dimensionedScalar n_;
        const dimensionedScalar rb_;
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

        
    	volScalarField f(); //hindering fuunnnnctiiion
    	dimensionedVector vStokes(); //Terminal velocity model
    	
    	volTensorField UHdiff();


public:

    //- Runtime type information
    TypeName("simpleP");


    // Constructors

        //- Construct from components
        simpleP
        (
            const dictionary& dict,
            const incompressibleTwoPhaseInteractingMixture& mixture,
            const word& modelName = ""
        );


    //- Destructor
    ~simpleP();


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
