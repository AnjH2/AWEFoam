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
    Foam::capillaryPressureModels::noPc

Description
     noPc J function for capillary pressure
     large portions of the model implementation is taken from openFuelCell2 
SourceFiles
    noPc.C

\*---------------------------------------------------------------------------*/

#ifndef noPc_H
#define noPc_H

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
                           Class noPc Declaration
\*---------------------------------------------------------------------------*/

class noPc
:
    public capillaryPressureModel
{
protected:

    // Protected data

        //- Dictionary
        dictionary noPcCoeffsSub_;

        


public:

    //- Runtime type information
    TypeName("noPc");


    // Constructors

        //- Construct from components
        noPc
        (
            const volScalarField& alphaWetting,
            const dictionary& dict,
            const word modelName=typeName
        );


    //- Destructor
    ~noPc()
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
