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
    Foam::speciesMixture

Description
    A two-phase species mixture model

SourceFiles
    speciesMixture.C

\*---------------------------------------------------------------------------*/

#ifndef speciesMixture_H
#define speciesMixture_H

#include "volFields.H"
#include "physicoChemicalConstants.H"
//#include "../porousProperties/porousProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class speciesMixture Declaration
\*---------------------------------------------------------------------------*/

class speciesMixture
{
protected:

    // Protected data
    	/*autoPtr<porousProperties> poroMPtr_ 
    	(
    		porousProperties
		(
		mesh
		)
    	);*/
	
	dictionary	species2Coeffs_;
	dictionary	species1Coeffs_;
	List<word> species2;
	List<word> species1;
	
	PtrList <dimensionedScalar> z_;
	PtrList <dimensionedScalar> Ea2_;
	PtrList <dimensionedScalar> D2_ref_;
	PtrList <dimensionedScalar> T2_ref_;
	PtrList <dimensionedScalar> MW_;
	PtrList <volScalarField> C2_;
	PtrList <volScalarField> D2_;
	PtrList <volScalarField> D2_eff_;
	//PtrList <volScalarField> D1_eff_;
	
	PtrList <volScalarField> C1_;
	//
	//start defining mixture coefficients
	//charge -constant
	//Mobility -constant
	dimensionedScalar	u_OH_;
	dimensionedScalar	u_K_;
	//phase 1
	dimensionedScalar		D1_C_star_;
	dimensionedScalar		D1_a_;
	dimensionedScalar		D1_b_;
	dimensionedScalar		D1_zeta_;
	dimensionedScalar		D1_ref_;
	//volScalarField			Delta_x1_;
	
	
        
        volScalarField	T_;
	
	//declearing external fields -> i belive i dont have to.
	//const porousProperties& poroM_(poroMPtr_());
	const volScalarField& epsilon_;
	const volScalarField& tau_;
	//const dimensionedScalar& R_;
	volScalarField			dummy_; //need a dummy who is equal 0.5
	volScalarField			D1_eff_;
	//making an object from class
	
public:

    // Constructors

        //- Construct from components
        speciesMixture
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const word& phase1NamE_,
            const word& phase2NamE_
        );


    //- Destructor
    //~speciesMixture()
    //{}
	virtual ~speciesMixture() = default;

    // Member Functions

        //- Return const-access to the eltroChemical class
       /* const porousProperties& poroM() const
        {
            return poroM_;
        }*/

        //making access to constant coefficients._______________________________
        //charge
        const PtrList <dimensionedScalar>& z() const
        {
            return z_;
        }
        

        //mobility
        const dimensionedScalar& u_OH() const
        {
            return u_OH_;
        }
        const dimensionedScalar& z_K() const
        {
            return u_K_;
        }
        //molarweight
        const PtrList <dimensionedScalar>& MW() const
        {
            return MW_;
        }


        //make access for fields__________________________________
        const volScalarField& T() const
        {
            return T_;
        }
        volScalarField& T()
        {
            return T_;
        }

        //- Return the species distribution of phase 2
        const PtrList <volScalarField>& C2() const
        {
            return C2_;
        }
        PtrList <volScalarField>& C2()
        {
            return C2_;
        }
        
        //- Return the species distribution of phase 1
        const PtrList <volScalarField>& C1() const
        {
            return C1_;
        }
        PtrList <volScalarField>& C1()
        {
            return C1_;
        }
        
        const PtrList <volScalarField>& D2_eff() const
        {
            return D2_eff_;
        }
        
        const volScalarField& D1_eff() const
        {
            return D1_eff_;
        }
       virtual tmp<volScalarField> DeltaX1( volScalarField X1) const
        {
            return D1_zeta_*pow(6*D1_C_star_-5,2)*(D1_a_*X1)/(1+D1_b_*X1);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
