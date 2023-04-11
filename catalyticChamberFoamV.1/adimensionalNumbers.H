//mean gas viscosity

volScalarField mum
(
    IOobject
    (
        "mum",
        runTime.timeName(),
        mesh
    ),
    (e1*tau1*thermo1.mu() + e2*tau2*thermo2.mu())/(et*taut)
);

//mean gas density

volScalarField rhom = (e1*thermo1.rho() + e2*thermo2.rho())/et;

// Reynold numbers: 
// Re0=U1*rhom/mum*as    catalyst  bed
// Re1=U1*dpa*rhom/mum   particles

volScalarField Re1
(
    IOobject
    (
        "Re1",
        runTime.timeName(),
        mesh
    ),
    e1*mag(U1)*rhom/(mum*as)
);

volScalarField Re2
(
    IOobject
    (
        "Re2",
        runTime.timeName(),
        mesh
    ),
    (e1*mag(U1)*dpa*rhom)/mum
);

//alpha of solid matrix

volScalarField DTc
(
    IOobject
    (
        "DTc",
        runTime.timeName(),
        mesh
    ),
    kC/(rhoC*cpC)
);

// mean gas rhomCpmf

Info<< "Solving rhomCpmf\n" << endl;

volScalarField rhomCpmf
(
    IOobject
    (
        "rhomCpmf",
        runTime.timeName(),
        mesh
    ),
    (e1*thermo1.rho()*thermo1.Cp() + e2*thermo2.rho()*thermo2.Cp())/et
);

// mean gas kappam 

Info<< "Solving fluid kappam\n" << endl;

volScalarField kappam
(
    IOobject
    (
        "kappam",
        runTime.timeName(),
        mesh
    ),
    (e1*thermo1.kappa() + e2*thermo2.kappa())/et
);

// Prandtl number
//Pr=mu/alphaeff=rho*cp/rho*kappa 

Info<< "Solving Pr\n" << endl;

volScalarField Pr
(
    IOobject
    (
        "Pr",
        runTime.timeName(),
        mesh
    ),
    (mum*rhomCpmf)/(rhom*kappam)
);

//Nulsselt number
//Nu=1.625Re^1/2Pr^1/3  RanzMarshall

Info<< "Solving Nu\n" << endl;

volScalarField Rep = Foam::pow(Re2, 0.6);

volScalarField Nu
(
    IOobject
    (
        "Nu",
        runTime.timeName(),
        mesh
    ),
    scalar(2) + 1.1*Rep*cbrt(Pr)
);

Info<< "Solving hh\n" << endl;

//Heat Transfer coefficient
//h=ktNu/dpa

volScalarField hh
(
    IOobject
    (
        "hh",
        runTime.timeName(),
        mesh
    ),
    kappam*Nu/dpa
);

Info << "adimensional Numbers: " << "Re1(max, min): " << max(Re1).value() << "," << min(Re1).value() << "/" << "Re2(max, min): " << max(Re2).value() << "," << min(Re2).value() << "/"<< "Pr(max,min): " << max(Pr).value() << "," << min(Pr).value() << "/" << "Nu(max, min):" << max(Nu).value() << "," << min(Nu).value() << "." << endl;

