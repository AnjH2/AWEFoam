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
            "porousProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    porousProperties_
    (
        IOobject
        (
            "porousProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    electrodes({"Ne","Pe"}),
    asL_(electrodes.size()),
    D_pore_(electrodes.size()),
    KL_(electrodes.size()+1),
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
       PeC_
    (
    	IOobject
        (
            "PeC",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless,Zero)
    ),
    
    NeC_
    (
    	IOobject
        (
            "NeC",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless,Zero)
    ),
        as_
    (
        IOobject
        (
            "as",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "as",
            dimless/dimLength,
            0
        )
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
            IOobject::NO_WRITE
        ),
        (Pe_+Ne_+VSMALL)*sigma_s_ref_*pow(1-epsilon_,tau_)//supress solid conductivity in all regions but electrode regions.
    ),
    K_
    (
    	IOobject
        (
            "K",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("K", dimArea,1e6)
    ),
    Kr_
    (
    	IOobject
        (
            "Kr",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("Kr", dimless,1)
    )
    
{
forAll(electrodes,i)
	{
		asL_.set
    	(
        	i,
        	new dimensionedScalar("as_"+electrodes[i], dimensionSet (  0, -1, 0, 0, 0, 0, 0),*this)
        );
        D_pore_.set
    	(
        	i,
        	new dimensionedScalar("D_pore_"+electrodes[i], dimensionSet (  0, 1, 0, 0, 0, 0, 0),*this)
        );
        KL_.set
    	(
        	i,
        	new dimensionedScalar("K_"+electrodes[i], dimensionSet (  0, 2, 0, 0, 0, 0, 0),*this)
        );
	}
        KL_.set
    	(
        	2,
        	new dimensionedScalar("K_Mem", dimensionSet (  0, 2, 0, 0, 0, 0, 0),*this)
        );
    //K_=Pe_*KL_[1]+Ne_*KL_[0]+Mem_*KL_[2]+(PeC_+NeC_)*dimensionedScalar(dimensionSet(  0, 2, 0, 0, 0, 0, 0),1e6);
    
    
        if (!K_.headerOk())
    {
    K_=Pe_*KL_[1]+Ne_*KL_[0]+Mem_*KL_[2]+(PeC_+NeC_)*dimensionedScalar(dimensionSet(  0, 2, 0, 0, 0, 0, 0),1e6);
    K_.correctBoundaryConditions();
        
    }

    if (!as_.headerOk())
    {
    
    as_ =
        dimensionedScalar
        (
            "as_Ne",
            dimless/dimLength,
            porousProperties_.get<scalar>("as_Ne")
        )*Ne_
      + dimensionedScalar
        (
            "as_Pe",
            dimless/dimLength,
            porousProperties_.get<scalar>("as_Pe")
        )*Pe_
      + dimensionedScalar("as_min", dimless/dimLength, SMALL);
        
    }
        if (!Kr_.headerOk())
    {
    
   Kr_ =
        dimensionedScalar
        (
            "Kr_Ne",
            dimless,
            porousProperties_.get<scalar>("Kr_Ne")
        )*Ne_
      + dimensionedScalar
        (
            "Kr_Pe",
            dimless,
            porousProperties_.get<scalar>("Kr_Pe")
        )*Pe_
        + dimensionedScalar
        (
            "Kr_Mem",
            dimless,
            porousProperties_.get<scalar>("Kr_Mem")
        )*Mem_
      + dimensionedScalar("Kr_min", dimless, 1)*(1-Ne_-Pe_-Mem_);
        
    }


}
// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //




// ************************************************************************* //
