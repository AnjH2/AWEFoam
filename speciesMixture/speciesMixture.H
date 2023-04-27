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

        //word phase1NamE_;
        //word phase2NamE_;

        volScalarField C_H2_1_;
        volScalarField C_O2_1_;
        volScalarField C_OH_1_;
        volScalarField C_H2O_1_;
        
        
        volScalarField C_H2_2_;
        volScalarField C_O2_2_;


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
