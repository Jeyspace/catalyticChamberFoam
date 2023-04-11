/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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
    catalyticChamberFoam

Description
    Solver for catalytic chambers of monopropellant thrusters with Two chemical reactions , one for a surface gas (heterogeneous reaction ) and another for interstitial gas ( homogeneous reaction ), mass Transfers for porous media with mass difusion coeficient calculate of Chagmant - Esnokg Theory . Mass and Heat transport relations is based on Wakao. For the doctoral theses: 
"ROBLEDO ASENCIO, J. C. Estudo numérico do processo de decomposição em câmaras catalíticas de propulsores monopropelentes. 2020. 189 p. IBI: <8JMKD3MGP3W34R/3UGHFGH>. (sid.inpe.br/mtc-m21c/2019/12.04.15.16-TDI). Tese (Doutorado em Combustão e Propulsão) - Instituto Nacional de Pesquisas Espaciais (INPE), São José dos Campos, 2019. Disponível em: <ibi:8JMKD3MGP3W34R/3UGHFGH>." and previously work: Asencio, J. R., Savonov, R. I., & Marques, R. I.. (2020). An Open-source Solver to Model the Catalytic Decomposition of Monopropellants for Space Thrusters. Journal of Aerospace Technology and Management, 12(J. Aerosp. Technol. Manag., 2020 12), e1120. https://doi.org/10.5028/jatm.v12.1111.



\*---------------------------------------------------------------------------*/

#include "fvCFD.H"			//reactingFoam
#include "turbulentFluidThermoModel.H"	//reactingFoam
#include "psiCombustionModel.H"		//reactingFoam
#include "psiReactionThermo.H"		//chemFoam
#include "psiChemistryModel.H"		//chemFoam
#include "chemistrySolver.H"		//chemFoam
#include "OFstream.H"			//chemFoam
#include "thermoPhysicsTypes.H"		//chemFoam
#include "basicMultiComponentMixture.H"	//chemFoam
#include "multivariateScheme.H"		//reactingFoam
#include "pimpleControl.H"		//reactingFoam
#include "pressureControl.H"		//reactingFoam
#include "fvOptions.H"			//reactingFoam
#include "localEulerDdtScheme.H"	//reactingFoam
#include "fvcSmooth.H"			//reactingFoam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"		  //reactingFoam/chemFoam

    #include "setRootCase.H"		  //reactingFoam/chemFoam
    #include "createTime.H"		  //reactingFoam/chemFoam
    #include "createMesh.H"		  //reactingFoam
    #include "createControl.H"		  //reactingFoam
    #include "createTimeControls.H"	  //reactingFoam
    #include "initContinuityErrs.H"	  //reactingFoam
    #include "createFields.H"		  //reactingFoam
    #include "createFieldRefs.H"	  //reactingFoam
    #include "createTransportProperties.H"//Transport Properties for Diffusion
    #include "createFvOptions.H"	  //reactingFoam
    #include "diffusionBinary.H"	  //Binary Diffusion reference

    turbulence->validate();

    if (!LTS)
    {
        #include "compressibleCourantNo.H" //reactingFoam
        #include "setInitialDeltaT.H"	   //reactingFoam
    }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H" //reactingFoam

        if (LTS)
        {
            #include "setRDeltaT.H"  //reactingFoam
        }
        else
        {
            #include "compressibleCourantNo.H"	//reactingFoam

            #include "setDeltaT.H"	        //reactingFoam/chemFoam

        }

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

            #include "adimensionalNumbers.H"
	    #include "diffusionCoefficient1.H"	
	    #include "diffusionCoefficient2.H"
            #include "massTransfer1.H"
	    #include "HeatTransfer1.H"
            #include "massTransfer2.H"
            #include "HeatTransfer2.H"
            #include "HeatTransferc.H"
            #include "HeatTransferM.H"
            #include "rho1Eqn.H"	//reactingFoam
            #include "rho2Eqn.H"	//reactingFoam

        while (pimple.loop())
        {
            #include "U1Eqn.H"		//reactingFoam
            #include "Y1Eqn.H"		//reactingFoam
            #include "E1Eqn.H"		//reactingFoam
 

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                if (pimple.consistent())
                {
                    #include "pc1Eqn.H"	//reactingFoam
                }
                else
                {
                    #include "p1Eqn.H"	//reactingFoam
                }

            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
	    
            #include "solveChemistry.H"	//chemFoam
            #include "Y2Eqn.H"		//chemFoam
            #include "E2Eqn.H"		//chemFoam
            #include "TcEqn.H"		//chemFoam

            #include "output.H"		//chemFoam

        }

        rho1 = thermo1.rho();
        rho2 = thermo2.rho();

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
