/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2016 OpenFOAM Foundation
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
    Foam::concentrationLossModel

Description

SourceFiles
    concentrationLossModel.C

\*---------------------------------------------------------------------------*/

#ifndef concentrationLossModel_H
#define concentrationLossModel_H

#include "fvCFD.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class concentrationLossModel Declaration
\*---------------------------------------------------------------------------*/

class concentrationLossModel
{
    // Private Member Functions

        //- No copy construct
        concentrationLossModel(const concentrationLossModel&) = delete;

        //- No copy assignment
        void operator=(const concentrationLossModel&) = delete;


protected:

    // Protected data
	const dictionary concentrationLossModelDict_;
	List<word> species2;
	PtrList <volScalarField> CR_Ne_;
	PtrList <volScalarField> CR_Pe_;

public:

    //- Runtime type information
    TypeName("concentrationLossModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        concentrationLossModel,
        dictionary,
        (const dictionary& dict,
         const fvMesh& mesh),
        (dict, mesh)
    );


    // Constructors

        //- Construct from components
        concentrationLossModel
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    // Selector
    static autoPtr<concentrationLossModel> New
    (
        const dictionary& dict,
        const fvMesh& mesh
    );


    //- Destructor
    virtual ~concentrationLossModel();


    // Member Functions
    	const PtrList <volScalarField> CR_Ne() const 
    	{
    		return CR_Ne_;
    	}
    	const PtrList <volScalarField> CR_Pe() const 
    	{
    		return CR_Pe_;
    	}
        //- Update the concentrationLoss
        virtual void correct(const PtrList <volScalarField>& C2_s_)=0;
        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
