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
    Foam::capillaryMassCentered

Description
    capillaryMassCentered relative velocity model,
    it defines properties for each electrode.
    lift only

SourceFiles
    capillaryMassCentered.C

\*---------------------------------------------------------------------------*/

#ifndef capillaryMassCentered_H
#define capillaryMassCentered_H

#include "relativeVelocityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace relativeVelocityModels
{

/*---------------------------------------------------------------------------*\
                         Class capillaryMassCentered Declaration
\*---------------------------------------------------------------------------*/

class capillaryMassCentered
:
    public relativeVelocityModel
{
	const incompressibleTwoPhaseInteractingMixture& mixture_;
	List<word> electrodes;
	const dictionary& dict_;
	
	volVectorField Urel_;
	
	
	const meshObjects::gravity& g_;
	const volScalarField& eps_;
	const volScalarField& K_;
	const volScalarField& kr_;
	const volScalarField& p_;
	const volScalarField& Solid_;
	
    // Private functions??
    volScalarField Mc();
    volScalarField Md();
    volScalarField Mm();
    volScalarField rhoS();
    volScalarField wd();
    volScalarField wc();
    
    



public:

    //- Runtime type information
    TypeName("capillaryMassCentered");


    // Constructors

        //- Construct from components
        capillaryMassCentered
        (
            const dictionary& dict,
            const incompressibleTwoPhaseInteractingMixture& mixture,
            const word& modelName = ""
        );


    //- Destructor
    ~capillaryMassCentered();


    // Member Functions
        //- Update the relative velocity
        void correctRelativeVelocity();
        //- Update the diffusion velocity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace relativeVelocityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
