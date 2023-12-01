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
#include "IOdictionary.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class porousProperties Declaration
\*---------------------------------------------------------------------------*/

class porousProperties
:
    public IOdictionary
{
protected:

    // Protected data
    	List<word> electrodes;	
    	PtrList <dimensionedScalar> as_;//electrode surface area.
    	PtrList <dimensionedScalar> D_pore_;//poreDiameter
	const volScalarField epsilon_;
	const volScalarField tau_;


	

	//regions of porouse media.
	const volScalarField Pe_;
	const volScalarField Ne_;
	const volScalarField Mem_;

	
	const dimensionedScalar sigma_s_ref_;
	const volScalarField sigma_s_eff_;


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
        const volScalarField& sigma_s_eff() const
        {
        	return sigma_s_eff_;
        }
        const volScalarField& Pe() const
        {
            return Pe_;
        }
        
        const volScalarField& Ne() const
        {
            return Ne_;
        }
        
        const PtrList <dimensionedScalar>& as() const
        {
            return as_;
        }
        
        const volScalarField& Mem() const
        {
            return Mem_;
        }
        const PtrList <dimensionedScalar>& D_pore() const
        {
            return D_pore_;
        }



        

        
        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
