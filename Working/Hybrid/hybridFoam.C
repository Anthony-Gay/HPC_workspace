/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

Application
    hybridFoam 

Description
    Call dsmc and rhoCentralFoam solvers in succsession

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "regionProperties.H"
#include "dsmcCloudMod2.H"
#include "psiThermo.H"
#include "turbulentFluidThermoModel.H"
#include "directionInterpolateMod.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #define NO_CONTROL
    //#define CREATE_MESH createMeshesPostProcess.H
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createTimeControls.H"

    #include "createMeshes.H"
    Info<< "Meshes created" << nl;

    // Create Fluid Fields
    #include "createFluidFields.H"
    #include "createFluidFieldsRefs.H"
    turbulence->validate();
    #include "readFluxSchemeMod.H"

    dimensionedScalar v_zero("v_zero", dimVolume/dimTime, 0.0);
    
    // Initialise particle regions
    #include "initialiseDSMC.H"

    // Creating Particle Fields
    #include "createDSMCFields.H"



    // Courant numbers used to adjust the time-step
    scalar CoNum = 0.0;
    scalar meanCoNum = 0.0;
    
     while (runTime.loop())
    {
       //objectRegistry().lookupObject()
       Info<< "Time = " << runTime.timeName() << nl << endl;
        Info<< nl << "SOLVING PARTICLE DOMAIN" << endl;
        // Solve Particle
        dsmcSolve.evolve();
        dsmcSolve.info();
        runTime.write();  

        // Solve Fluid
        Info<< nl << "SOLVING FLUID DOMAIN" << endl;
        #include "solveFluid.H"

          

        // Output runtime
        dimensionedScalar simTime=runTime.time();
         dimensionedScalar simDt=runTime.deltaT();
         dimensionedScalar endTime=runTime.endTime();
         scalar remainIts=(endTime.value()/simDt.value())-(simTime.value()/simDt.value());
         scalar timePerIt=runTime.elapsedCpuTime()/(simTime.value()/simDt.value());

        Info<< nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << "  Iterations: " << (simTime.value()/simDt.value()) << "/" << (endTime.value()/simDt.value())
            << "  Projected time to finish: " << (remainIts*timePerIt)/60 << " min"
            << endl;       
             
             
    }

    return 0;
}


// ************************************************************************* //
