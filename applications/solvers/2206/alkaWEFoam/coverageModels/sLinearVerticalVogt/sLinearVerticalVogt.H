/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2015 OpenFOAM Foundation
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
    Foam::sLinearVerticalVogt

Description
    sLinearVerticalVogt coverage model

SourceFiles
    sLinearVerticalVogt.C

\*---------------------------------------------------------------------------*/

#ifndef sLinearVerticalVogt_H
#define sLinearVerticalVogt_H

#include "volFields.H"
#include "unitConversion.H"
#include "coverageModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace coverageModels
{

/*---------------------------------------------------------------------------*\
                         Class linear Declaration
\*---------------------------------------------------------------------------*/

class sLinearVerticalVogt
:
    public coverageModel
{
    // Private data

        //- exchange current density -> read
        List<word> electrodes;
        const dictionary& dict_;
        const volScalarField& J_;
        const volScalarField& T_;
        const dimensionedScalar J_lim_;
        const dimensionedScalar J_scale_;
        const dimensionedScalar T_ref_;
        PtrList <dimensionedScalar> CD_;    //drag coefficient for bubbles
        PtrList <dimensionedScalar> c_; //linearScale that converts coverage to saturation
	    PtrList <dimensionedScalar> rd_;
	    const PtrList <dimensionedScalar>& as_;
        
        
        const volScalarField& Pe_;
        const volScalarField& Ne_; 
        const volVectorField& U_;
        const dimensionedScalar& rhoc_;
        const volScalarField& rhod_;
        const meshObjects::gravity& g_;
        //- volume fraction of gas -> read
        const volScalarField& alphad_;
        
        // functions
        volScalarField correction();
        
        
        
public:

    //- Runtime type information
    TypeName("sLinearVerticalVogt");


    // Constructors

        //- Construct from components
        sLinearVerticalVogt
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const incompressibleTwoPhaseInteractingMixture& mixture,
        const word& modelName = ""
        );


    //- Destructor
    ~sLinearVerticalVogt();


    // Member Functions

        //- Update coverage
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace coverageModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
