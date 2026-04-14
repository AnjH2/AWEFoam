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
    under the terms of the GNU generalR Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU generalR Public License
    for more details.

    You should have received a copy of the GNU generalR Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::sherwoodModels::generalR

Description
     sherwood correction model for generic newtonian flow.

SourceFiles
    generalR.C

\*---------------------------------------------------------------------------*/

#ifndef generalR_H
#define generalR_H

#include "sherwoodModel.H"
#include "relativeVelocityModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace sherwoodModels
{

/*---------------------------------------------------------------------------*\
                           Class generalR Declaration
\*---------------------------------------------------------------------------*/

class generalR
:
    public sherwoodModel
{
protected:

    // Protected data

        //- Dictionary
        dictionary generalRCoeffsSub1_;
        //- Maximum sherwood
        //dimensionedScalar muMax_;
	PtrList <dimensionedScalar> d_;
	PtrList <dimensionedScalar> abc_;

	



public:

    //- Runtime type information
    TypeName("generalR");


    // Constructors

        //- Construct from components
        generalR
        (
            const word& name,
            const incompressibleTwoPhaseInteractingMixture& mixture,
            const dictionary& sherwoodPropertiesSub1,
            const word modelName=typeName
        );


    //- Destructor
    ~generalR()
    {}


    // Member Functions
	const PtrList <dimensionedScalar>& d()const 
	{
		return d_;
	}
	tmp<volScalarField> ki(const int i);

        //- Read transportProperties dictionary
        bool read(const dictionary& sherwoodPropertiesSub1);
        

	const tmp<volScalarField> sh(const int i)
	{
		return 1+abc_[0]*pow(ReR(d_[0])*(Ne_+NeC_)+ReR(d_[1])*(Pe_+PeC_)+ReR((d_[1]+d_[0])/2)*Mem_,abc_[1])*pow(Sc(i),abc_[2]);
	}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace sherwoodModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
