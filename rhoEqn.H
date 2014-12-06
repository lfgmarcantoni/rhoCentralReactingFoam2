  // surfaceScalarField phi("phi", aphiv_pos*rho_pos + aphiv_neg*rho_neg);
phi = aphiv_pos*rho_pos + aphiv_neg*rho_neg;

Info<<"rho prev max/min"<<max(rho)<<""<<min(rho)<<endl;
//Info<<"rho prev max/min"<<max(rho)<<""<<min(rho)<<endl;
//volScalarField rhon("rhon",rho.dimensions(),rho.value()); 
volScalarField rhon = rho; 
     // solving continuity 
        solve(fvm::ddt(rhon) + fvc::div(phi));
rho=rhon; 
rhon.correctBoundaryConditions(); 
//fvOptions.correct(rho); 

Info<<"rho post max/min"<<max(rho)<<""<<min(rho)<<endl;

//thermo.correct();
rho=rhon; 
