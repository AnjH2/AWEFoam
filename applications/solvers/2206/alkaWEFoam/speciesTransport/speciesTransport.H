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
    Foam::speciesTransport

Description
    A two-phase species mixture model

SourceFiles
    speciesTransport.C

\*---------------------------------------------------------------------------*/

#ifndef speciesTransport_H
#define speciesTransport_H


#include "volFields.H"
#include "physicoChemicalConstants.H"
#include "massAndSpeciesTransferModel.H"

//#include "IOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class speciesTransport Declaration
\*---------------------------------------------------------------------------*/
class viscosityModel;
class speciesTransport
//:
   //public IOdictionary
{
protected:

    // Protected data
    	const massAndSpeciesTransferModel& mSTaPtr_;
	dictionary	speciesTransportCoeffs_;
	//wire radius, is either constructed from booleans or read from 0/rb. for the sherwoodModel its converted into diameter
    volScalarField wr_;
	autoPtr<sherwoodModel> shModelW_;
	List<word> species2;
	PtrList <volScalarField> C2_s_;
	const PtrList <volScalarField>& C2_;
	const volScalarField& alphal_;
	const volScalarField& Pe_;
	const volScalarField& Ne_;
	const volScalarField& PeC_;
	const volScalarField& NeC_;
	const volScalarField& Mem_;
	const volScalarField& as_;
	const PtrList <dimensionedScalar>& MW_;
	
	

	





public:
	//TypeName("speciesTransport");
    // Constructors

        //- Construct from components
        speciesTransport
        (
    		const fvMesh& mesh,
    		const dictionary& dict,
    		const massAndSpeciesTransferModel& mSTaPtr
        );


    //- Destructor
 	virtual ~speciesTransport() = default;


    // Member Functions
    const PtrList <volScalarField>& C2_s() const
        {
            return C2_s_;
        }
    
	void correct(const int i, const volScalarField& theta_);
}; 


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
