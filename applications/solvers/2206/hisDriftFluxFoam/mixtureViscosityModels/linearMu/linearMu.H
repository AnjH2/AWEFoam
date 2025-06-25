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
    Foam::mixtureViscosityModels::linearMu

Description
     Viscosity correction model for generic newtonian flow.

SourceFiles
    linearMu.C

\*---------------------------------------------------------------------------*/

#ifndef linearMu_H
#define linearMu_H

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
                           Class linearMu Declaration
\*---------------------------------------------------------------------------*/

class linearMu
:
    public mixtureViscosityModel
{
protected:

    // Protected data

        //- Dictionary
        dictionary linearMuCoeffsSub1_;
        dictionary linearMuCoeffsSub2_;
        //- Maximum viscosity
        //dimensionedScalar muMax_;

	PtrList <dimensionedScalar> mud_ref_;
	PtrList <dimensionedScalar> T_muref_;
	PtrList <dimensionedScalar> n_mu_;
	PtrList <volScalarField> mud_;
	volScalarField a_;
	volScalarField mud_m_;
	        
	        
        //- linearMu phase fraction
        const volScalarField& alpha_;
        //const volScalarField& rhod_;


public:

    //- Runtime type information
    TypeName("linearMu");


    // Constructors

        //- Construct from components
        linearMu
        (
            const word& name,
            const dictionary& viscosityPropertiesSub1,
            const dictionary& viscosityPropertiesSub2,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const word modelName=typeName
        );


    //- Destructor
    ~linearMu()
    {}


    // Member Functions

        //- Return the mixture viscosity
        //  given the viscosity of the continuous phase
        tmp<volScalarField> mu(const volScalarField& rhod) const;
	const PtrList<volScalarField>& mud() const
	{
		return mud_;
	}
	
	virtual	const volScalarField& mud_m() const
	{
		return mud_m_;
	}
	void mud_m_correct(const volScalarField& p_water,const dimensionedScalar& p_num);

        //- Read transportProperties dictionary
        bool read(const dictionary& viscosityPropertiesSub1);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace mixtureViscosityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
