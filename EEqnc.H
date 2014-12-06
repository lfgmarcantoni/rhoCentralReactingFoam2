{


//Info<<"heName"<<he.name()<<endl ; 

    volScalarField hec =he ;

Info<<"heName"<<he.name()<<endl;

    fvScalarMatrix EEqnc
    (
        fvm::ddt(rho, hec) + fvm::div(phi, hec)
      + fvc::ddt(rho, K) + fvc::div(phi, K)
//      + (
//            hec.name() == "e"
//          ? fvc::div
 //           (
 //               fvc::absolute(phi/fvc::interpolate(rho), U),
 //               p,
  //              "div(phiv,p)"
   //         )
    //      : -dpdt
      //  )
//      - fvm::laplacian(turbulence->alphaEff(), he)
//     ==
//        reaction->Sh()
 //     + fvOptions(rho, hec)
    );

    EEqnc.relax();

    fvOptions.constrain(EEqnc);

    EEqnc.solve();

    fvOptions.correct(hec);

    thermo.correct();

    Info<< "min/max(T) after correction= "
        << min(T).value() << ", " << max(T).value() << endl;
}
