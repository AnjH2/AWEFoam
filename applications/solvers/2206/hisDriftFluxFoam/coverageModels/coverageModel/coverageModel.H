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
    Foam::coverageModel

Description

SourceFiles
    coverageModel.C

\*---------------------------------------------------------------------------*/

#ifndef coverageModel_H
#define coverageModel_H

#include "fvCFD.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class coverageModel Declaration
\*---------------------------------------------------------------------------*/

class coverageModel
{
    // Private Member Functions

        //- No copy construct
        coverageModel(const coverageModel&) = delete;

        //- No copy assignment
        void operator=(const coverageModel&) = delete;


protected:

    // Protected data
	const dictionary coverageModelDict_;
	const dictionary transportPropertiesDict_;
	volScalarField	theta_;
	dimensionedScalar p_num_;

public:

    //- Runtime type information
    TypeName("coverageModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        coverageModel,
        dictionary,
        (const dictionary& dict,
         const fvMesh& mesh),
        (dict, mesh)
    );


    // Constructors

        //- Construct from components
        coverageModel
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    // Selector
    static autoPtr<coverageModel> New
    (
        const dictionary& dict,
        const fvMesh& mesh
    );


    //- Destructor
    virtual ~coverageModel();


    // Member Functions


        //- Return the diffusion velocity of the dispersed phase
        const volScalarField& theta() const
        {
            return theta_;
        }
        
        void uP_num(const dimensionedScalar p_num) {
        	p_num_=p_num;
        }

        //- Update the coverage
        virtual void correct() = 0;
        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
