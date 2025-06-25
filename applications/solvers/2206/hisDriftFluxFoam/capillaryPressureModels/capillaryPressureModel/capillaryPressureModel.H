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
    Foam::capillaryPressureModels

Description
    A namespace for incompressible capillaryPressureModel implementations.

Class
    Foam::capillaryPressureModel

Description
    An abstract base class for incompressible capillaryPressureModels.

SourceFiles
    capillaryPressureModel.C
    capillaryPressureModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef capillaryPressureModel_H
#define capillaryPressureModel_H

#include "fvCFD.H"
#include "dictionary.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "dimensionedScalar.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class capillaryPressureModel Declaration
\*---------------------------------------------------------------------------*/

class capillaryPressureModel
{

protected:

    // Protected data
        const volScalarField& alphaWetting_;
        const dictionary& dict_;
        const volScalarField& Solid_;
        
        volScalarField pc_;
        volScalarField dpcds_;

    // Private Member Functions

        //- No copy construct
        capillaryPressureModel(const capillaryPressureModel&) = delete;

        //- No copy assignment
        void operator=(const capillaryPressureModel&) = delete;


public:

    //- Runtime type information
    TypeName("capillaryPressureModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            capillaryPressureModel,
            dictionary,
            (
                const volScalarField& alphaWetting,
                const dictionary& dict
            ),
            (alphaWetting,dict)
        );


    // Selectors

        //- Return a reference to the selected viscosity model
        static autoPtr<capillaryPressureModel> New
        (
            const volScalarField& alphaWetting,
            const dictionary& dict
        );


    // Constructors

        //- Construct from components
        capillaryPressureModel
        (
            const volScalarField& alphaWetting,
            const dictionary& dict
        );


    //- Destructor
    virtual ~capillaryPressureModel() = default;



    // Member Functions
    
    //- return reference functions
    const volScalarField& pc() const
    {
        return pc_;
    }

    const volScalarField& dpcds() const
    {
        return dpcds_;
    }
    //- Correct the capillary preussure and its derivative
    virtual void correct() = 0;

    //- Solve for saturation
    //virtual void updateSaturation() = 0; dont know if i need it ANJ



        //- Read transportProperties dictionary
        virtual bool read(const dictionary& viscosityPropertiesSub1) = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
