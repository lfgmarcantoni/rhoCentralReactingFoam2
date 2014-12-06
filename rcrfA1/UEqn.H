//discretizacion de div(rho u u) según kurganov

/*
 surfaceVectorField rhoU_pos =
            fvc::interpolate(rhoU, pos, "reconstruct(U)");
        surfaceVectorField rhoU_neg =
            fvc::interpolate(rhoU, neg, "reconstruct(U)");
*/
Info<<"solving  Ueqn"<<endl; 
 surfaceVectorField phiUp 
(
            (aphiv_pos*rhoU_pos + aphiv_neg*rhoU_neg)
          + (a_pos*p_pos + a_neg*p_neg)*mesh.Sf()
);
//    fvVectorMatrix UEqn
  solve  (
        fvm::ddt(rhoU)
      + fvc::div(phiUp)
    );
//fvOptions.correct(U);
//K=0.5*magSqr(U);  // se actualiza K 

//    UEqn.relax();
//    UEqn.solve(); 
U = rhoU/rho;
//  U.dimensionedInternalField() =
//            rhoU.dimensionedInternalField()
//           /rho.dimensionedInternalField();
 U.correctBoundaryConditions();
  rhoU.boundaryField() = rho.boundaryField()*U.boundaryField();
  rhoU.dimensionedInternalField() = rho.dimensionedInternalField()*U.dimensionedInternalField();
 volScalarField rhoBydt(rho/runTime.deltaT());

fvOptions.correct(U); 
K=0.5*magSqr(U);  // se actualiza K 

if (thermo.dpdt())
{
//    dpdt = fvc::DDt(p);

//DpDt = fvc::DDt(surfaceScalarField("phiUp", phi/fvc::interpolate(rho)), p);

}


