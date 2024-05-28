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
    Foam::sherwoodModels::general

Description
     sherwood correction model for generic newtonian flow.

SourceFiles
    general.C

\*---------------------------------------------------------------------------*/

#ifndef general_H
#define general_H

#include "sherwoodModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace sherwoodModels
{

/*---------------------------------------------------------------------------*\
                           Class general Declaration
\*---------------------------------------------------------------------------*/

class general
:
    public sherwoodModel
{
protected:

    // Protected data

        //- Dictionary
        dictionary generalCoeffsSub1_;
        //- Maximum sherwood
        //dimensionedScalar muMax_;
	const dimensionedScalar d_;
	PtrList <dimensionedScalar> abc_;

	



public:

    //- Runtime type information
    TypeName("general");


    // Constructors

        //- Construct from components
        general
        (
            const word& name,
            const incompressibleTwoPhaseInteractingMixture& mixture,
            const dictionary& sherwoodPropertiesSub1,
            const word modelName=typeName
        );


    //- Destructor
    ~general()
    {}


    // Member Functions
	const dimensionedScalar& d()const {
		return d_;
	}
	tmp<volScalarField> ki(const int i);

        //- Read transportProperties dictionary
        bool read(const dictionary& sherwoodPropertiesSub1);
        

	const tmp<volScalarField> sh(const int i){
		return 1+abc_[0]*pow(Re(d_),abc_[1])*pow(Sc(i),abc_[2]);
	}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace sherwoodModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
