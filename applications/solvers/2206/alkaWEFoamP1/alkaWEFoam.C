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
#include "sherwoodModel.H"
#include "turbulenceModel.H"
#include "CompressibleTurbulenceModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "gaussLaplacianScheme.H"
#include "uncorrectedSnGrad.H"
#include "IOobjectList.H"
#include "loopControl.H"
#include "condKOH.H"
#include "m_KOH.H"

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
    const volScalarField& rho1 = mixture.rhod();
    PtrList <volScalarField>& C2 = speciesP.C2();
    const PtrList <dimensionedScalar>& MW = speciesP.MW();
    const PtrList <dimensionedScalar>& z = speciesP.z();
    const volScalarField& T = speciesP.T();

    
    relativeVelocityModel& UdmModel(UdmModelPtr());
    
    coverageModel& thetaModel(thetaModelPtr());
    const volScalarField& theta=thetaModel.theta();
    
    concentrationLossModel& CLModel(CLModelPtr());
    const volScalarField& CRa_Ne = CLModel.CRa_Ne();
    const volScalarField& CRa_Pe = CLModel.CRa_Pe();
    const volScalarField& CRc_Ne = CLModel.CRc_Ne();
    const volScalarField& CRc_Pe = CLModel.CRc_Pe();
    
    massAndSpeciesTransferModel& mSTa(mSTaPtr());
    
    const volScalarField& epsilon1 = poroM.eps();
    const PtrList <dimensionedScalar>& as = poroM.as();
    const volScalarField& Ne = poroM.Ne();
    const volScalarField& Pe = poroM.Pe();
    const volScalarField& Mem = poroM.Mem();
    
    const dimensionedScalar& F=Foam::constant::physicoChemical::F;
    const dimensionedScalar& R=Foam::constant::physicoChemical::R;
    
   
    //listing species for species transport equations.

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;
    
    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "driftCourantNo.H"
        #include "setDeltaT.H"

        ++runTime;
        
        
        Info<< "Time = " << runTime.timeName() << nl << endl;
        mSTa.correct_waterPartialPressure(cal_mKOH(T,C2[3]/1000)*dimensionedScalar(dimensionSet(1,0,0,0,-1,0,0),1));
        
        // --- initialising mixture velocity
       
        mixture.correct(mSTa.p_water());
        
        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
           // mSTa.correct_waterPartialPressure(cal_mKOH(T,C2[3]/1000)*dimensionedScalar(dimensionSet(1,0,0,0,-1,0,0),1));	
            #include "alphaControls.H"
            
            UdmModel.correct();
            
            thetaModel.uP_num(mixture.p_num());
            thetaModel.correct();
            
            #include "alphaEqnSubCycle.H"
	    //mixture.correct_rhoc();

	    mixture.correct_rhod(mSTa.p_water());

            mixture.correct(mSTa.p_water());

            //mixture.correct_D1();

            mixture.correct_nu();

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
            
            speciesP.correct_DOH((calKappa(T,C2[3]/1000)));
            
            speciesP.correct_D2_eff(1-alpha2);
            
            react.correct(mSTa.p_water(),mSTa.p_water_pure(),mixture.p_num(),cal_mKOH(T,C2[3]/1000));
            
            for (int k=0; k<=outerChemicalCorrections; k++)
            {
            
            	//speciesP.correct_DOH((calKappa(T,C2[3]/1000)));
            	
            	//speciesP.correct_D2_eff(1-alpha2);
            	
            	CLModel.correct(sTp.C2_s());
            	
            	for (int j=0; j<=potentialCorrections; j++)
            	{
            	
            		#include "potentialEqn.H"
            	}

            	for (int j=0; j<=speciesCorrections; j++)
            	{
            		#include "CiEqn.H"
            		//mixture.correct();
            	}
	    }
        }
	if (runTime.writeTime())
	{
		
        	INe = poroM.sigma_s_eff()*fvc::grad(UNe);
        	INe.correctBoundaryConditions();
        	IPe = poroM.sigma_s_eff()*fvc::grad(UPe);
        	INe.correctBoundaryConditions();
	}
        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
