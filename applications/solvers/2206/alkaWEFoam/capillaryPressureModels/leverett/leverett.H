/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014 OpenFOAM Foundation
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
    Foam::capillaryPressureModels::leverett

Description
     Leverett J function for capillary pressure
     large portions of the model implementation is taken from openFuelCell2 
SourceFiles
    leverett.C

\*---------------------------------------------------------------------------*/

#ifndef leverett_H
#define leverett_H

#include "capillaryPressureModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class incompressibleTwoPhaseInteractingMixture;

namespace capillaryPressureModels
{

/*---------------------------------------------------------------------------*\
                           Class leverett Declaration
\*---------------------------------------------------------------------------*/

class leverett
:
    public capillaryPressureModel
{
protected:

    // Protected data

        //- Dictionary
        dictionary leverettCoeffsSub_;
        //- Surface tension
        dimensionedScalar sigma_;

        //- Contact angle
        scalar theta_;

        // irreducible saturation
        dimensionedScalar alphaWetMin_;

        // maximal saturation
        dimensionedScalar alphaWetMax_;

        //- effective saturation of phase 2
        volScalarField alphaWetEff_;
        
        
        
        //- porosity
        const volScalarField& eps_;
        const volScalarField& K_;
        const volScalarField& Solid_;
        


public:

    //- Runtime type information
    TypeName("leverett");


    // Constructors

        //- Construct from components
        leverett
        (
            const volScalarField& alphaWetting,
            const dictionary& dict,
            const word modelName=typeName
        );


    //- Destructor
    ~leverett()
    {}


    // Member Functions
        //updates pc
        virtual void correct();
        
        
        //- Read transportProperties dictionary
        bool read(const dictionary& viscosityPropertiesSub1);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace capillaryPressureModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
