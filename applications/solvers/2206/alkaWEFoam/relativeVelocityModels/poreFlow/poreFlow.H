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
    Foam::poreFlow

Description
    poreFlow relative velocity model,
    it defines properties for each electrode.
    diffusion and lift only

SourceFiles
    poreFlow.C

\*---------------------------------------------------------------------------*/

#ifndef poreFlow_H
#define poreFlow_H

#include "relativeVelocityModel.H"
#include "coverageModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace relativeVelocityModels
{

/*---------------------------------------------------------------------------*\
                         Class poreFlow Declaration
\*---------------------------------------------------------------------------*/

class poreFlow
:
    public relativeVelocityModel
{
	autoPtr<relativeVelocityModel> continuousModel_;
	autoPtr<relativeVelocityModel> dispersedModel_;
	List<word> electrodes;
	volTensorField Ddm_;
	const volScalarField& alphaEff_;//needs a capillaryPressure model! !!! and the wetting phase is always alphac !!!
	
	
	volScalarField chi(); //smooth function to blend capillarity and dispersed flow.
    



public:

    //- Runtime type information
    TypeName("poreFlow");


    // Constructors

        //- Construct from components
        poreFlow
        (
            const dictionary& dict,
            const incompressibleTwoPhaseInteractingMixture& mixture,
            const word& modelName = ""
        );


    //- Destructor
    ~poreFlow();


    // Member Functions

        //- Update the diffusion velocity
        virtual void correct();
        
        	virtual const volTensorField& Ddm() const
        {
            return Ddm_;
        }
        virtual bool runDiffusion()
        {
        	return dispersedModel_->runDiffusion() or continuousModel_->runDiffusion();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace relativeVelocityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
