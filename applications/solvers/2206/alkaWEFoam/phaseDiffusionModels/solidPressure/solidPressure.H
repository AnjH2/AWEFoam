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
    Foam::solidPressure

Description
    Porous-region bubble-dispersion closure with an anisotropic baseline
    coefficient and an exponential gas-saturation augmentation. The latter
    mimics the solid-pressure-like increase in dispersion at high gas
    saturation used in the associated electrolyzer model.

SourceFiles
    solidPressure.C


\*---------------------------------------------------------------------------*/

#ifndef solidPressure_H
#define solidPressure_H

#include "phaseDiffusionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseDiffusionModels
{

/*---------------------------------------------------------------------------*\
                         Class solidPressure Declaration
\*---------------------------------------------------------------------------*/

class solidPressure
:
    public phaseDiffusionModel
{
    // Private data
    
        List<word> electrodes;
        const dictionary& dict_;
        const meshObjects::gravity& g_;
        dimensionedScalar Dp_;
        dimensionedScalar Dn_;
        //const dimensionedScalar n_;
        PtrList <dimensionedScalar> rd_;
        vector ep_; //unit vector parallel to g
        volTensorField D_;
	
	    const scalar  ANe_;
	    const scalar  CbNe_;
	    const scalar  APe_;
	    const scalar  CbPe_;
        const scalar  alpha0_;

        
    	volScalarField f(); //hindering fuunnnnctiiion
    	volVectorField vStokes(); //Terminal velocity model
    	
    	volTensorField UHdiff();


public:

    //- Runtime type information
    TypeName("solidPressure");


    // Constructors

        //- Construct from components
        solidPressure
        (
            const dictionary& dict,
            const incompressibleTwoPhaseInteractingMixture& mixture,
            const word& modelName = ""
        );


    //- Destructor
    ~solidPressure();


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
