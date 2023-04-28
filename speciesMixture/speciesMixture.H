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

	dictionary	speciesCoeffs_;

	//start defining mixture coefficients
	//charge -constant
	dimensionedScalar	z_OH_;
	dimensionedScalar	z_H2O_;
	dimensionedScalar	z_H2_;
	dimensionedScalar	z_O2_;
	//Mobility -constant
	dimensionedScalar	u_OH_;
	dimensionedScalar	u_K_;
	
	//diffusion related constants -constant
	dimensionedScalar	Ea_H2O_;
	dimensionedScalar	Ea_OH_;
	dimensionedScalar	Ea_H2_;
	dimensionedScalar	Ea_O2_;
	
	dimensionedScalar	D_H2O_ref_;
	dimensionedScalar	D_OH_ref_;
	dimensionedScalar	D_H2_ref_;
	dimensionedScalar	D_O2_ref_;
		//constant
	dimensionedScalar	T_ref_H2_;
	dimensionedScalar	T_ref_O2_;
	dimensionedScalar	T_ref_H2O_;
	dimensionedScalar	T_ref_OH_;
	
	//mole weight -constant
	dimensionedScalar	MW_H2_;
	dimensionedScalar	MW_O2_;
	dimensionedScalar	MW_OH_;
	dimensionedScalar	MW_H2O_;
	dimensionedScalar	MW_K_;
	
	//declearing fields. Phase1
        volScalarField C_H2_1_;
        volScalarField C_O2_1_;
        volScalarField C_OH_1_;
        volScalarField C_H2O_1_;
        
        //phase2
        volScalarField C_H2_2_;
        volScalarField C_O2_2_;
	
	//declearing external fields -> i belive i dont have to.
	const volScalarField& T_;
	const volScalarField& epsilon_;
	const volScalarField& tau_;
	//const dimensionedScalar& R_;
	
	
	//derived fields
	volScalarField	D_H2O;
	volScalarField	D_OH;
	volScalarField	D_H2;
	volScalarField	D_O2;
	
	volScalarField	D_H2O_eff;
	volScalarField	D_OH_eff;
	volScalarField	D_H2_eff;
	volScalarField	D_O2_eff;


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
    ~speciesMixture()
    {}


    // Member Functions

        /*const word& phase1NamE() const
        {
            return phase1NamE_;
        }

        const word& phase2NamE() const
        {
            return phase2NamE_;
        }*/

        //making access to constant coefficients._______________________________
        //charge
        const dimensionedScalar& z_OH() const
        {
            return z_OH_;
        }
        const dimensionedScalar& z_H2O() const
        {
            return z_H2O_;
        }
        const dimensionedScalar& z_H2() const
        {
            return z_H2_;
        }
        const dimensionedScalar& z_O2() const
        {
            return z_O2_;
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
        const dimensionedScalar& MW_OH() const
        {
            return MW_OH_;
        }
        const dimensionedScalar& MW_H2O() const
        {
            return MW_H2O_;
        }
        const dimensionedScalar& MW_H2() const
        {
            return MW_H2_;
        }
        const dimensionedScalar& MW_O2() const
        {
            return MW_O2_;
        }
        const dimensionedScalar& MW_K() const
        {
            return MW_K_;
        }
        //make access for fields__________________________________

        //- Return the phase-fraction of phase 1
        const volScalarField& C_H2_1() const
        {
            return C_H2_1_;
        }

        //- Return the phase-fraction of phase 1
        volScalarField& C_H2_1()
        {
            return C_H2_1_;
        }
	
	//- Return the phase-fraction of phase 1
        const volScalarField& C_O2_1() const
        {
            return C_O2_1_;
        }

        //- Return the phase-fraction of phase 1
        volScalarField& C_O2_1()
        {
            return C_O2_1_;
        }
        
        //- Return the phase-fraction of phase 1
        const volScalarField& C_OH_1() const
        {
            return C_OH_1_;
        }

        //- Return the phase-fraction of phase 1
        volScalarField& C_OH_1()
        {
            return C_OH_1_;
        }
        
        //- Return the phase-fraction of phase 1
        const volScalarField& C_H2O_1() const
        {
            return C_H2O_1_;
        }

        //- Return the phase-fraction of phase 1
        volScalarField& C_H2O_1()
        {
            return C_H2O_1_;
        }
        
        //- Return the phase-fraction of phase 2
        const volScalarField& C_H2_2() const
        {
            return C_H2_2_;
        }

        //- Return the phase-fraction of phase 2
        volScalarField& C_H2_2()
        {
            return C_H2_2_;
        }
        
        //- Return the phase-fraction of phase 2
        const volScalarField& C_O2_2() const
        {
            return C_O2_2_;
        }

        //- Return the phase-fraction of phase 2
        volScalarField& C_O2_2()
        {
            return C_O2_2_;
        }
        

        
        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
