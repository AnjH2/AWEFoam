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
    along with OpenFOAM.  If not, see <http://(www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "porousProperties.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::porousProperties::porousProperties
(
    const fvMesh& mesh
)
:
    IOdictionary
    (
        IOobject
        (
            "transportProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    electrodes({"Ne","Pe"}),
    as_(electrodes.size()),
    epsilon_
    (
        IOobject
        (
            "eps",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    ),
    
    tau_
    (
        IOobject
        (
            "tau",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    ),
    Pe_
    (
    	IOobject
        (
            "Pe",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    
    Ne_
    (
    	IOobject
        (
            "Ne",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    
    Mem_
    (
    	IOobject
        (
            "Mem",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
   
        sigma_s_ref_
	(
		"sigma_s_ref",
		dimensionSet ( -1, -3, 3, 0, 0, 2,0),
		*this
	),
    sigma_s_eff_
    (
    	IOobject
        (
            "sigma_s_eff",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        (Pe_+Ne_+Mem_*1e-6)*sigma_s_ref_*pow(1-epsilon_,tau_)//supress solid conductivity in all regions but electrode regions.
    )
{
forAll(electrodes,i)
	{
		as_.set
    	(
        	i,
        	new dimensionedScalar("as_"+electrodes[i], dimensionSet (  0, -1, 0, 0, 0, 0, 0),*this)
        );
	}

}
// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //




// ************************************************************************* //
