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
    Foam::speciesProperties

Description
    Species-property container for the alkaline electrolyte model. Reads
    concentrations and molecular properties, evaluates temperature-dependent
    diffusivities, and applies porous/gas corrections to effective diffusion.

SourceFiles
    speciesProperties.C

\*---------------------------------------------------------------------------*/

#ifndef speciesProperties_H
#define speciesProperties_H

#include "fvCFD.H"
#include "volFields.H"
#include "physicoChemicalConstants.H"
#include "IOdictionary.H"
//#include "../porousProperties/porousProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class speciesProperties Declaration
\*---------------------------------------------------------------------------*/

class speciesProperties
:
    public IOdictionary
{
protected:

	word phase2NamE_;
	dictionary	species2Coeffs_;
	List<word> species2;
	
	PtrList <dimensionedScalar> z_;
	PtrList <dimensionedScalar> Ea2_;
	PtrList <dimensionedScalar> D2_ref_;
	PtrList <dimensionedScalar> T2_ref_;
	PtrList <dimensionedScalar> MW_;
	PtrList <volScalarField> C2_;
	PtrList <volScalarField> D2_;
	PtrList <volScalarField> D2_eff_;
	PtrList <volScalarField> t_;
	volScalarField  D2_ambi_;
	volScalarField  C2_T_;
	volScalarField  C2_0_;
	dimensionedScalar	u_OH_;
	dimensionedScalar	u_K_;
	volScalarField	T_;
	const volScalarField& epsilon_;
	const volScalarField& tau_;
	const volScalarField& Mem_;
	
public:

    // Constructors

        //- Construct from components
        speciesProperties
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    //~speciesProperties()
    //{}
	virtual ~speciesProperties() = default;

    // Member Functions


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
        const dimensionedScalar& u_K() const
        {
            return u_K_;
        }
        
        //transference numbers
        const PtrList <volScalarField>& t() const
        {
            return t_;
        }
        //conductivity
        
        
        //molarweight
        const PtrList <dimensionedScalar>& MW() const
        {
            return MW_;
        }
        const volScalarField kappa() const
        {
            return -1*pow(Foam::constant::physicoChemical::F,2)/(Foam::constant::physicoChemical::R*T_)*(z_[4]*D2_[4]-z_[3]*D2_[3])*z_[3]*C2_[3];
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
        const volScalarField& C2_T() const
        {
            return C2_T_;
        }
        const volScalarField& C2_0() const
        {
            return C2_0_;
        }        
        const PtrList <volScalarField>& D2_eff() const
        {
            return D2_eff_;
        }
        const PtrList <volScalarField>& D2() const
        {
            return D2_;
        }
        const volScalarField& D2_ambi() const
        {
            return D2_ambi_;
        }
        void correct_DOH(volScalarField kappa)
        {

        	D2_[3]=kappa*Foam::constant::physicoChemical::R*T_/(2*C2_[3]*pow(Foam::constant::physicoChemical::F*z_[3],2));
        }
	void correct_D2_eff(volScalarField alphad_)
        {
        	forAll(species2,i)
        		{
        			D2_eff_[i]=Foam::pow(epsilon_,tau_)*Foam::pow(1-min(alphad_,0.99)*(1-Mem_),1.5)*D2_[i];
        		}


            D2_ambi_=D2_eff_[3];
        }
     void correct();   

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
