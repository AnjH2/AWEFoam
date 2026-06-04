/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2015 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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

Namespace
    Foam::sherwoodModels

Description
    A namespace for incompressible sherwoodModel implementations.

Class
    Foam::sherwoodModel

Description
    An abstract base class for incompressible sherwoodModels.

SourceFiles
    sherwoodModel.C
    sherwoodModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef sherwoodModel_H
#define sherwoodModel_H


#include "fvCFD.H"
#include "dictionary.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "dimensionedScalar.H"
#include "runTimeSelectionTables.H"
#include "incompressibleTwoPhaseInteractingMixture.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class sherwoodModel Declaration
\*---------------------------------------------------------------------------*/

class sherwoodModel

{

protected:

    // Protected data

        word name_;
        const incompressibleTwoPhaseInteractingMixture& mixture_;
        List<word> species2;
        List<word> electrodes;
        dictionary sherwoodPropertiesSub1_;
        
        const PtrList <volScalarField>& D2_;
        const volVectorField& Udm_;

	const volScalarField& Pe_;
	const volScalarField& PeC_;
	const volScalarField& Ne_;
	const volScalarField& NeC_;
	const volScalarField& Mem_;
	const volScalarField& lCh_;
    // Private Member Functions

        //- No copy construct
        sherwoodModel(const sherwoodModel&) = delete;

        //- No copy assignment
        void operator=(const sherwoodModel&) = delete;


public:

    //- Runtime type information
    TypeName("sherwoodModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            sherwoodModel,
            dictionary,
            (
                const word& name,
                const incompressibleTwoPhaseInteractingMixture& mixture,
                const dictionary& sherwoodPropertiesSub1,
                const volScalarField& lCh //characteristic length
            ),
            (name, mixture, sherwoodPropertiesSub1, lCh)
        );


    // Selectors

        //- Return a reference to the selected viscosity model
        static autoPtr<sherwoodModel> New
        (
            const word& name,
            const incompressibleTwoPhaseInteractingMixture& mixture,
            const dictionary& sherwoodPropertiesSub1,
            const volScalarField& lCh
        );


    // Constructors

        //- Construct from components
        sherwoodModel
        (
            const word& name,
            const incompressibleTwoPhaseInteractingMixture& mixture,
            const dictionary& sherwoodPropertiesSub1,
            const volScalarField& lCh
        );


    //- Destructor
    virtual ~sherwoodModel() = default;

    // Member Functions
    
	virtual  scalar s() const = 0;
	virtual tmp<volScalarField> ki(const int i) = 0;
    
    
	tmp<volScalarField> Re(const volScalarField d)const{
		return (mag(mixture_.U())*d/mixture_.nuc());
	}
	tmp<volScalarField> ReR(const volScalarField d)const{
		return (mag(mixture_.rho()/(mixture_.rhoc()*mixture_.alpha2())*Udm_)*d/mixture_.nuc());
	}
	
	tmp<volScalarField> Sc(const int i)const{
		return (mixture_.nuc()/D2_[i]);
	}
	virtual const tmp<volScalarField> sh(const int i) = 0;
	
        //- Read transportProperties dictionary
        virtual bool read(const dictionary& sherwoodPropertiesSub1) = 0;

	
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
