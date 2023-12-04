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

Application
    driftFluxFoam

Group
    grpMultiphaseSolvers

Description
    Solver for two incompressible fluids using the mixture approach with
    the drift-flux approximation for relative motion of the phases.

    Used for simulating the settling of the dispersed phase and other
    similar separation problems.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "CMULES.H"
#include "subCycle.H"

#include "porousProperties.H"
#include "speciesProperties.H"
#include "incompressibleTwoPhaseInteractingMixture.H"
#include "reactionProperties.H"
#include "relativeVelocityModel.H"
#include "coverageModel.H"
#include "concentrationLossModel.H"
#include "speciesTransport.H"
#include "massAndSpeciesTransferModel.H"
#include "turbulenceModel.H"
#include "CompressibleTurbulenceModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "gaussLaplacianScheme.H"
#include "uncorrectedSnGrad.H"
#include "IOobjectList.H"
#include "loopControl.H"
#include "condKOH.H"
#include "Tensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Solver for two incompressible fluids using the mixture approach with"
        " the drift-flux approximation for relative motion of the phases.\n"
        "Used for simulating the settling of the dispersed phase and other"
        " similar separation problems."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createTimeControls.H"
    #include "createFields.H"
    #include "physicoChemicalConstants.H"

    
    
    #include "initContinuityErrs.H"

    volScalarField& alpha2(mixture.alpha2());
    const dimensionedScalar& rho2 = mixture.rhoc();
    PtrList <volScalarField>& C2 = speciesP.C2();
    PtrList <volScalarField>& C1 = speciesP.C1();
    const PtrList <dimensionedScalar>& MW = speciesP.MW();
    const PtrList <dimensionedScalar>& z = speciesP.z();
    const volScalarField& T = speciesP.T();

    
    relativeVelocityModel& UdmModel(UdmModelPtr());
    
    coverageModel& thetaModel(thetaModelPtr());
    const volScalarField& theta=thetaModel.theta();
    
    concentrationLossModel& CLModel(CLModelPtr());
    const PtrList <volScalarField>& CR_Ne = CLModel.CR_Ne();
    const PtrList <volScalarField>& CR_Pe = CLModel.CR_Pe();
    
    massAndSpeciesTransferModel& mSTa(mSTaPtr());
    
    const volScalarField& epsilon1 = poroM.eps();
    const PtrList <dimensionedScalar>& as = poroM.as();
    const volScalarField& Ne = poroM.Ne();
    const volScalarField& Pe = poroM.Pe();
    //const volScalarField& Mem = poroM.Mem();
    
    const dimensionedScalar& F=Foam::constant::physicoChemical::F;
    const dimensionedScalar& R=Foam::constant::physicoChemical::R;
    
   
    //listing species for species transport equations.

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;
    mSTa.correct_waterPartialPressure();
    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;
        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "alphaControls.H"
            UdmModel.correct();
            //thetaModel.correct();
            thetaModel.correct();
            #include "alphaEqnSubCycle.H"
            
	    //mixture.correct_rhoc();
	    mixture.correct_rhod();
            mixture.correct();
            //mixture.correct_D1();
            mixture.correct_nu();
            speciesP.correct_D2_eff(1-alpha2);
            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
		
                #include "pEqn.H"
            }
            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
            CLModel.correct(sTp.C2_s());

            for (int j=0; j<=potentialCorrections; j++)
            {
            #include "potentialEqn.H"
            }

            for (int j=0; j<=speciesCorrections; j++)
            {
            #include "CiEqn.H"
            }

        }
	if (runTime.writeTime())
	{
		
        	INe = poroM.sigma_s_eff()*fvc::grad(UNe);
        	IPe = poroM.sigma_s_eff()*fvc::grad(UPe);

        	Ie = -(calKappa(T,C2[3]/1000)*(pow(epsilon1,poroM.tau())*pow(alpha2,1.5))*fvc::grad(Ue));
        	
        }
        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
