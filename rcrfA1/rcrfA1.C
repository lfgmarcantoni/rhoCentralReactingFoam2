/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    reactingFoam

Description
    Solver for combustion with chemical reactions.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "turbulenceModel.H"
#include "zeroGradientFvPatchFields.H"
#include "fixedRhoFvPatchScalarField.H"
#include "psiCombustionModel.H"
#include "multivariateScheme.H"
//#include "pimpleControl.H"
#include "fvIOoptionList.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readGravitationalAcceleration.H"
    #include "createFields1.H"
    #include "createFvOptions.H"
//    #include "initContinuityErrs.H"
    #include "readTimeControls.H"
//    scalar CoNum = 0.0005; 
    #include "readFluxScheme.H"

//    #include "compressibleCourantNo.H"
  //  #include "setInitialDeltaT.H"

    //pimpleControl pimple(mesh);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


//#   include "readFluxScheme.H"

    dimensionedScalar v_zero("v_zero",dimVolume/dimTime, 0.0);

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
//#       include "readTimeControls.H"
//#       include "setInitialDeltaT.H"
//#       include "setDeltaT.H"


  //      #include "setDeltaT.H"

//volScalarField Gamma(thermo.gamma()); 
//thermo.correct();
// thermo.correct();
Info<<"gamma max/min" <<max(Gamma)<<min(Gamma)<<endl;
        surfaceScalarField rho_pos
	(
            fvc::interpolate(rho, pos, "reconstruct(rho)")
	);
            surfaceScalarField rho_neg
        (
            "rho_neg",
            fvc::interpolate(rho, neg, "reconstruct(rho)")
        );

        surfaceVectorField U_pos
        (
            "U_pos",
            fvc::interpolate(U, pos, "reconstruct(U)")
        );

        surfaceVectorField U_neg
        (
            "U_neg",
            fvc::interpolate(U, neg, "reconstruct(U)")
        );

Info<<"se calculo Upos y Uneg"<<endl; 

        volScalarField rPsi (1.0/psi);
   surfaceScalarField rPsi_pos
        (
            "rPsi_pos",
            fvc::interpolate(rPsi, pos, "reconstruct(T)")
        );
       
     surfaceScalarField rPsi_neg
        (
            "rPsi_neg",
            fvc::interpolate(rPsi, neg, "reconstruct(T)")
        );
 
//        surfaceVectorField U_pos("U_pos", rhoU_pos/rho_pos);
  //      surfaceVectorField U_neg("U_neg", rhoU_neg/rho_neg);


 surfaceVectorField rhoU_pos("rhoU_pos", U_pos*rho_pos);
 surfaceVectorField rhoU_neg("rhoU_neg", U_neg*rho_neg);


//     surfaceScalarField Gamma_pos =
//            fvc::interpolate(Gamma, pos, "reconstruct(T)");
//        surfaceScalarField Gamma_neg =
//            fvc::interpolate(Gamma, neg, "reconstruct(T)");

        surfaceScalarField p_pos("p_pos", rho_pos*rPsi_pos);
        surfaceScalarField p_neg("p_neg", rho_neg*rPsi_neg);
        surfaceScalarField phiv_pos("phiv_pos", U_pos & mesh.Sf());
        surfaceScalarField phiv_neg("phiv_neg", U_neg & mesh.Sf());
 

//volScalarField  R=1.0/(psi*thermo.T());

// volScalarField c = sqrt(thermo.Cp()/thermo.Cp()*rPsi);
volScalarField Cp = thermo.Cp(); 
volScalarField Cv = thermo.Cv(); 

c = sqrt(Gamma*rPsi);  //sqrt((Cp/Cv)*(Cp-Cv)*T);
volScalarField c2(sqrt(thermo.Cp()/thermo.Cv()*rPsi));

 surfaceScalarField cSf_pos
        (
            "cSf_pos",
            fvc::interpolate(c, pos, "reconstruct(T)")*mesh.magSf()
        );

        surfaceScalarField cSf_neg
        (
            "cSf_neg",
            fvc::interpolate(c, neg, "reconstruct(T)")*mesh.magSf()
        );

    surfaceScalarField ap
        (
            "ap",
            max(max(phiv_pos + cSf_pos, phiv_neg + cSf_neg), v_zero)
        );

       surfaceScalarField am
        (
            "am",
            min(min(phiv_pos - cSf_pos, phiv_neg - cSf_neg), v_zero)
        );

        surfaceScalarField a_pos("a_pos", ap/(ap - am));

	surfaceScalarField amaxSf("amaxSf", max(mag(am), mag(ap)));


       surfaceScalarField aSf("aSf", am*a_pos);

        if (fluxScheme == "Tadmor")
        {
            aSf = -0.5*amaxSf;
            a_pos = 0.5;
        }

        surfaceScalarField a_neg("a_neg",1.0 - a_pos);

        phiv_pos *= a_pos;
        phiv_neg *= a_neg;


        surfaceScalarField aphiv_pos("aphiv_pos",phiv_pos - aSf);
        surfaceScalarField aphiv_neg("aphiv_neg", phiv_neg + aSf);

        amaxSf = max(mag(aphiv_pos), mag(aphiv_neg));
//#       include "compressibleCourantNo.H"
        #include "compressibleCourantNo.H"
#       include "readTimeControls.H"
//#       include "setInitialDeltaT.H"
#       include "setDeltaT.H"


  //      #include "setDeltaT.H"

        runTime++;
        Info<< "Time = " << runTime.timeName() << nl << endl;
//rho = p*psi ; 
        #include "rhoEqn.H"
        #include "UEqn.H"
       #include "YEqn.H"
        if (he.name()=="e"){
        #include "EEqn2.H"
        }
        if (he.name()=="h") {
//       #include "HEqn.H"
Info<<"no yet implemented"<<endl; 
return 0; 
        }
//#include "rhoEqn.H"
//        #include "YEqn2.H"
//        #include"EEqnc.H"

//Gamma = thermo.gamma(); 
//#include "GammaEqn.H"
p.dimensionedInternalField() =
            rho.dimensionedInternalField()
           /psi.dimensionedInternalField();
        p.correctBoundaryConditions();
       rho.boundaryField() = psi.boundaryField()*p.boundaryField();
//thermo->
//        while (pimple.loop())
  //      {
         
            // --- Pressure corrector loop
    //        while (pimple.correct())
     //       {
         //       #include "pEqn.H"
       //     }

           // if (pimple.turbCorr())
           // {
//                turbulence->correct();
           // }
       //}

        runTime.write();
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
thermo.correct();
Gamma = thermo.gamma();  
M = mag(U)/c ; 
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
