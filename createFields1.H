Info<< "Creating reaction model\n" << endl;

autoPtr<combustionModels::psiCombustionModel> reaction
(
    combustionModels::psiCombustionModel::New(mesh)
);

psiReactionThermo& thermo = reaction->thermo();
thermo.validate(args.executable(), "h", "e");

basicMultiComponentMixture& composition = thermo.composition();
PtrList<volScalarField>& Y = composition.Y();

word inertSpecie(thermo.lookup("inertSpecie"));

Info <<"R gas "<<max(thermo.Cp()-thermo.Cv())<<"--"<<min(thermo.Cp()-thermo.Cv())<<endl;

volScalarField& p = thermo.p();
const volScalarField& psi = thermo.psi();
const volScalarField& T = thermo.T();
//volScalarField& Gamma = thermo.gamma();

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
p*psi
//    thermo.rho()
);



volScalarField R
(
    IOobject
    (
        "R",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
thermo.Cp()-thermo.Cv()
//    thermo.rho()
);


volScalarField Cp1
(
    IOobject
    (
        "Cp1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
thermo.Cp()
//    thermo.rho()
);



volScalarField Cv1
(
    IOobject
    (
        "Cv1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
thermo.Cv()
//    thermo.rho()
);



volScalarField T2
(
    IOobject
    (
        "T2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
thermo.T()
//    thermo.rho()

);

Info<<"density ini"<<max(rho)<<"--"<<min(rho)<<endl ; 
Info<<"Temp ini"<<max(T)<<"--"<<min(T)<<endl ; 
Info<<"press ini"<<max(p)<<"--"<<min(p)<<endl ; 
Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);



//volScalarField& p = thermo.p();
//const volScalarField& psi = thermo.psi();
//const volScalarField& T = thermo.T();
//volScalarField& Gamma = thermo.gamma();

surfaceScalarField pos
(
    IOobject
    (
        "pos",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("pos", dimless, 1.0)
);



surfaceScalarField neg
(
    IOobject
    (
        "neg",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("neg", dimless, -1.0)
);

volVectorField rhoU
(
    IOobject
    (
        "rhoU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*U
);


//#include "compressibleCreatePhi.H"
/*
Info << "Creating turbulence model.\n" << nl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);
*/
// Set the turbulence into the reaction model
//reaction->setTurbulence(turbulence());



Info<< "Creating field dpdt\n" << endl;
volScalarField dpdt
(
    IOobject
    (
        "dpdt",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("dpdt", p.dimensions()/dimTime, 0)
);

Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K("K", 0.5*magSqr(U));
volScalarField& he = thermo.he();
//volScalarField  rhohe = rho*he; 
//volScalarField rhoK = rho*K; 

volScalarField rhohe
(
    IOobject
    (
        "rhohe",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*he
);


volScalarField rhoE
(
    IOobject
    (
        "rhoE",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rhohe+rho*K
);


multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

forAll(Y, i)
{
    fields.add(Y[i]);
}
fields.add(thermo.he());

volScalarField dQ
(
    IOobject
    (
        "dQ",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("dQ", dimEnergy/dimTime, 0.0)
);

Info<<"Cp max/min"<<max(thermo.Cp())<<min(thermo.Cp())<<endl; 
Info<<"Cv max/min"<<max(thermo.Cv())<<min(thermo.Cv())<<endl; 

PtrList<volScalarField> rhoY(Y.size());

forAll (Y, i)
{
         rhoY.set
         (
             i,
             new volScalarField
             (
                 IOobject
                 (
                     "rhoY",
                     runTime.timeName(),
		     mesh,
                     IOobject::NO_READ,
                     IOobject::AUTO_WRITE
                 ),
                 rho*Y[i]
             )
         );
    	 fields.add(Y[i]);
Info<<"rhoY[i] "<<i<<rhoY[i]<<endl;//rhoY.name()<<endl; 
}

//PtrList<volScalarField> rhoY(Y*rho);


volScalarField c
(
    IOobject
    (
        "c",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("c", dimLength/dimTime, 0)
//    sqrt((thermo.Cp()/thermo.Cv())*(thermo.Cp()-thermo.Cv())*thermo.T())
);


volScalarField c2
(
    IOobject
    (
        "c2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    sqrt(thermo.Cp()/thermo.Cv()*(1/thermo.psi()))
);


volScalarField Gamma
(
    IOobject
    (
        "Gamma",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.gamma()
);


volScalarField rhoGamma
(
    IOobject
    (
        "rhoGamma",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*Gamma
);


volScalarField  energyVar
(
    IOobject
    (
        "energyVar",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    rhohe/rho
);

surfaceScalarField phi("phi", mesh.Sf() & fvc::interpolate(rhoU));


Info << "Creating turbulence model.\n" << nl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);

reaction->setTurbulence(turbulence());
/*
volScalarField DpDt
(
    IOobject
    (
        "DpDt",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("DpDt", p.dimensions()/dimTime, 0)
);
*/
//volScalarField DpDt =fvc::DDt(surfaceScalarField("phi", phi/fvc::interpolate(rho)), p);

volScalarField M
(
    IOobject
    (
        "M",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mag(U)/(sqrt(Gamma*(thermo.Cp()-thermo.Cv())*thermo.T()))
);

volScalarField  w = rhohe-p; 


//volScalarField hs = thermo.hs() ; 

volScalarField rhoT
(
    IOobject
    (
        "rhoT",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
//    mag(U)/(sqrt(Gamma*(thermo.Cp()-thermo.Cv())*thermo.T()))
rho*thermo.T()
);

//surfaceScalarField phiT("phi", mesh.Sf() & fvc::interpolate(rhoU));

surfaceScalarField  phiT
(
    IOobject
    (
        "phiT",
        runTime.timeName(),
        mesh
    ),
   mesh,
    dimensionedScalar("pos", dimVolume*dimTime/dimTemperature, 0.0)
);

