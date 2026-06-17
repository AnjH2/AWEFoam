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
    Foam::capillaryPressureModels::BrooksAndCorey

Description
     BrooksAndCorey J function for capillary pressure
     large portions of the model implementation is taken from openFuelCell2 
SourceFiles
    BrooksAndCorey.C

\*---------------------------------------------------------------------------*/

#ifndef BrooksAndCorey_H
#define BrooksAndCorey_H

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
                           Class BrooksAndCorey Declaration
\*---------------------------------------------------------------------------*/

class BrooksAndCorey
:
    public capillaryPressureModel
{
protected:

    // Protected data

        //- Dictionary
        dictionary BrooksAndCoreyCoeffsSub_;

        // entry capillary pressure
        volScalarField pc0_;

        // Exponent
        volScalarField beta_;
    
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
        const volScalarField& Mem_;
        


public:

    //- Runtime type information
    TypeName("BrooksAndCorey");


    // Constructors

        //- Construct from components
        BrooksAndCorey
        (
            const volScalarField& alphaWetting,
            const dictionary& dict,
            const word modelName=typeName
        );


    //- Destructor
    ~BrooksAndCorey()
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
