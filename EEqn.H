{
//    volScalarField& he = thermo.he();
//volScalarField  rhohe = rho*he; 
//volScalarField rhoK = rho*K; 

  surfaceScalarField he_pos 
(
"he_pos",
            fvc::interpolate(he, pos, "reconstruct(T)")
); 
  surfaceScalarField he_neg
(
"he_neg", 
         fvc::interpolate(he, neg, "reconstruct(T)")
); 
  surfaceScalarField K_pos(
"K_pos",
            fvc::interpolate(K, pos, "reconstruct(K)")
);
  surfaceScalarField K_neg(
"K_neg",
            fvc::interpolate(K, neg, "reconstruct(K)")
);


//surfaceScalarField rhoK_pos =  rho_pos*K_pos; 
//surfaceScalarField rhoK_neg =  rho_neg*K_neg;

//if (he.name() == "e"){
//Info<<"se resuelve la eq de energia total"<<endl; 
//}
surfaceScalarField phihe(
"phihe",
  aphiv_pos*(rho_pos*(he_pos)+p_pos)
+ aphiv_neg*(rho_neg*(he_neg)+p_neg)
+ aSf*p_pos-aSf*p_neg
); 
//surfaceScalarField phihe = aphiv_pos*((he_pos)+p_pos) + aphiv_neg*((he_neg)+p_neg)+aSf*p_pos-aSf*p_neg ;

// contiene las contribuciones de div(rhoEU) y  div(up) 
//Info<<"se resolvera la eq de la energia total"<<endl ; 


//if (he.name() == "h"){
//Info<<"se resuelve la ecuacion de la entalpia "<<endl; 
//phihe = phihe-aphiv_pos*p_pos-aphiv_neg*p_neg-aSf*p_pos+aSf*p_neg;
//}

//surfaceScalarField phiP =  phiv_pos*(p_pos) + phiv_neg*(p_neg);

surfaceScalarField phiK(
 "phiK",
  aphiv_pos*(rho_pos*(0.5*magSqr(U_pos))) + aphiv_neg*(rho_neg*(0.5*magSqr(U_neg)))
); 
//    fvScalarMatrix EEqn
 solve   (
        fvm::ddt(rhohe)  + fvc::div(phihe)
      + fvc::ddt(rhoK) + fvc::div(phiK)
//      - dpdt //(
 //         he.name() == "h"
//         ?  -dpdt
//          :  dpdtDummy
//      )
//      - fvm::laplacian(turbulence->alphaEff(), he)
  //  ==
   //     reaction->Sh()
//      + fvOptions(rhohe)
    );

//   EEqn.relax();

    //fvOptions.constrain(EEqn);

  //  EEqn.solve();

  //  fvOptions.correct(he);

//he.dimensionedInternalField() =
//rhohe.dimensionedInternalField()/rho.dimensionedInternalField();

he = rhohe/rho;
energyVar = rhohe /rho ; 
energyVar.boundaryField()= rhohe.boundaryField()/rho.boundaryField();

he.boundaryField()= rhohe.boundaryField()/rho.boundaryField();
he.correctBoundaryConditions();
energyVar.correctBoundaryConditions();

Info<<"max /  min "<<he.name()<<max(he)<<" " <<min(he)<<endl; 
//    EEqn.relax();

 //   fvOptions.constrain(EEqn);

 //   EEqn.solve();

    fvOptions.correct(he);

    thermo.correct();
//volScalarField Tt = thermo.T();
    Info<< "min/max(T) = "
        << min(T).value() << ", " << max(T).value() << endl;
}
