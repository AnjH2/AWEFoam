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
    Foam::mixtureViscosityModels::Beckermann

Description
     Viscosity correction model for generic newtonian flow.

SourceFiles
    Beckermann.C

\*---------------------------------------------------------------------------*/

#ifndef Beckermann_H
#define Beckermann_H

#include "mixtureViscosityModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class incompressibleTwoPhaseInteractingMixture;

namespace mixtureViscosityModels
{

/*---------------------------------------------------------------------------*\
                           Class Beckermann Declaration
\*---------------------------------------------------------------------------*/

class Beckermann
:
    public mixtureViscosityModel
{
protected:

    // Protected data

        //- Dictionary
        dictionary BeckermannCoeffsSub1_;
        dictionary BeckermannCoeffsSub2_;
        //- Maximum viscosity
        //dimensionedScalar muMax_;

	PtrList <dimensionedScalar> mud_ref_;
	PtrList <dimensionedScalar> T_muref_;
	PtrList <dimensionedScalar> n_mu_;
	PtrList <volScalarField> mud_;
	volScalarField a_;
	volScalarField mud_m_;
	dimensionedScalar rhoc_;
	dimensionedScalar kr_;
	        
	        
        //- Beckermann phase fraction
        const volScalarField& alpha_;
        //const volScalarField& rhod_;


public:

    //- Runtime type information
    TypeName("Beckermann");


    // Constructors

        //- Construct from components
        Beckermann
        (
            const word& name,
            const dictionary& viscosityPropertiesSub1,
            const dictionary& viscosityPropertiesSub2,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const word modelName=typeName
        );


    //- Destructor
    ~Beckermann()
    {}


    // Member Functions

        //- Return the mixture viscosity
        //  given the viscosity of the continuous phase
        tmp<volScalarField> mu(const volScalarField& muc, const volScalarField& rhod) const;

        //- Read transportProperties dictionary
        bool read(const dictionary& viscosityPropertiesSub1);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace mixtureViscosityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
