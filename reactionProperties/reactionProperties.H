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
    Foam::reactionProperties

Description
    A two-phase species mixture model

SourceFiles
    reactionProperties.C

\*---------------------------------------------------------------------------*/

#ifndef reactionProperties_H
#define reactionProperties_H


#include "volFields.H"
#include "physicoChemicalConstants.H"
#include "IOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class reactionProperties Declaration
\*---------------------------------------------------------------------------*/

class reactionProperties
:
    public IOdictionary
{
protected:

    // Protected data
	dictionary	reactionCoeffs_;
	List<word> electrodes;
	List<word> species2;
	PtrList <dimensionedScalar> C_ref_;
	PtrList <dimensionedScalar> C_sat_;
	PtrList <dimensionedScalar> psi_;
	PtrList <volScalarField> Psi_;
	PtrList <dimensionedScalar> E0_ref_;
	PtrList <dimensionedScalar> T0_E0_;
	PtrList <dimensionedScalar> Ds_;
	PtrList <volScalarField> E0_;
	PtrList <dimensionedScalar> i0_ref_;
	PtrList <dimensionedScalar> Ea_;
	PtrList <dimensionedScalar> T0_i0_;
	PtrList <volScalarField> i0_;
	PtrList <dimensionedScalar> alphaT_;
	const dimensionedScalar n_;	
	const dimensionedScalar u_K_;
	const dimensionedScalar u_OH_;
	dimensionedScalar t_OH_;
	dimensionedScalar tau_c_;
	const volScalarField&	T_;





public:
	//TypeName("reactionProperties");
    // Constructors

        //- Construct from components
        reactionProperties
        (
            const fvMesh& mesh
        );


    //- Destructor
 	virtual ~reactionProperties() = default;


    // Member Functions

        const PtrList <dimensionedScalar>& C2_ref() const
        {
            return C_ref_;
        }
         const PtrList <dimensionedScalar>& C_sat() const
        {
            return C_sat_;
        }

        const PtrList <dimensionedScalar>& psi() const
        {
            return psi_;
        }
	PtrList <volScalarField>& Psi()
        {
            return Psi_;
        }
        const PtrList <volScalarField>& i0() const
        {
            return i0_;
        }
        
        const PtrList <dimensionedScalar>& alphaT() const
        {
            return alphaT_;
        }
        
        const PtrList <volScalarField>& E0() const
        {
            return E0_;
        }
        const dimensionedScalar& t_OH() const
        {
            return t_OH_;
        }
        const dimensionedScalar& n() const
        {
            return n_;
        }
        
        const dimensionedScalar& tau_c() const
        {
            return tau_c_;
        }
      
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
