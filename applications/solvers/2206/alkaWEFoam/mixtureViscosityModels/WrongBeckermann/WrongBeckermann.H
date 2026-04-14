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
    Foam::mixtureViscosityModels::WrongBeckermann

Description
     Viscosity correction model for generic newtonian flow.

SourceFiles
    WrongBeckermann.C

\*---------------------------------------------------------------------------*/

#ifndef WrongBeckermann_H
#define WrongBeckermann_H

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
                           Class WrongBeckermann Declaration
\*---------------------------------------------------------------------------*/

class WrongBeckermann
:
    public mixtureViscosityModel
{
protected:

    // Protected data

        //- Dictionary
        dictionary WrongBeckermannCoeffsSub1_;
        dictionary WrongBeckermannCoeffsSub2_;
        //- Maximum viscosity
        //dimensionedScalar muMax_;


	dimensionedScalar mud_;
	dimensionedScalar rhoc_;
	        
	        
        //- WrongBeckermann phase fraction
        const volScalarField& alpha_;
        //const volScalarField& rhod_;


public:

    //- Runtime type information
    TypeName("WrongBeckermann");


    // Constructors

        //- Construct from components
        WrongBeckermann
        (
            const word& name,
            const dictionary& viscosityPropertiesSub1,
            const dictionary& viscosityPropertiesSub2,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const word modelName=typeName
        );


    //- Destructor
    ~WrongBeckermann()
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
