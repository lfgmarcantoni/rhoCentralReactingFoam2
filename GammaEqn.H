//  surfaceScalarField Gamma_pos =
//            fvc::interpolate(Gamma, pos, "reconstruct(he)");
 // surfaceScalarField Gamma_neg =
 //           fvc::interpolate(Gamma, neg, "reconstruct(he)");

   surfaceScalarField phiG("phiG", aphiv_pos*Gamma_pos*rho_pos + aphiv_neg*Gamma_neg*rho_neg);

//Info<<"rho prev max/min"<<max(rho)<<""<<min(rho)<<endl;
//Info<<"rho prev max/min"<<max(rho)<<""<<min(rho)<<endl;

     // solving continuity 
        solve(fvm::ddt(rhoGamma) + fvc::div(phiG));
Gamma=rhoGamma/rho; 
Gamma.correctBoundaryConditions(); 
fvOptions.correct(Gamma); 

Info<<"gamma post max/min"<<max(Gamma)<<""<<min(Gamma)<<endl;

