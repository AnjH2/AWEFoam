/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
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
    Foam::porousProperties

Description
    A two-phase species mixture model

SourceFiles
    porousProperties.C

\*---------------------------------------------------------------------------*/

#ifndef porousProperties_H
#define porousProperties_H


#include "volFields.H"
#include "physicoChemicalConstants.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class porousProperties Declaration
\*---------------------------------------------------------------------------*/

class porousProperties
{
protected:

    // Protected data

	const volScalarField epsilon_;
	const volScalarField tau_;


public:
	//TypeName("porousProperties");
    // Constructors

        //- Construct from components
        porousProperties
        (
            const fvMesh& mesh
        );


    //- Destructor
 	virtual ~porousProperties() = default;


    // Member Functions

        const volScalarField& eps() const
        {
            return epsilon_;
        }

        
        const volScalarField& tau() const
        {
            return tau_;
        }


        

        
        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
