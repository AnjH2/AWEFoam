/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

Global
    CourantNo

Description
    Calculates and outputs the mean and maximum Courant Numbers.

\*---------------------------------------------------------------------------*/

scalar CoNum = 0.0;
scalar meanCoNum = 0.0;

scalar DriftCoNum = 0.0;
scalar DriftmeanCoNum = 0.0;

{
    scalarField sumPhi
    (
        fvc::surfaceSum(mag(phi))().primitiveField()
    );

    CoNum = 0.5*gMax(sumPhi/mesh.V().field())*runTime.deltaTValue();

    meanCoNum =
        0.5*(gSum(sumPhi)/gSum(mesh.V().field()))*runTime.deltaTValue();
   
   scalarField sumUdm
    (
        fvc::surfaceSum(mag(phir))().primitiveField()
    );

    DriftCoNum = 0.5*gMax(sumUdm/mesh.V().field())*runTime.deltaTValue();

    DriftmeanCoNum =
        0.5*(gSum(sumUdm)/gSum(mesh.V().field()))*runTime.deltaTValue();

}

Info<< "Courant Number mean: " << meanCoNum
    << " max: " << CoNum << endl;
    
Info<< "Drift Courant Number mean: " << DriftmeanCoNum
    << " max: " << DriftCoNum << endl;   

// ************************************************************************* //
