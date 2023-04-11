IOdictionary porousProperties
(
    IOobject
    (
        "porousProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

scalar phi0 = readScalar(porousProperties.lookup("phis"));
scalar porosity2 = readScalar(porousProperties.lookup("porosityP"));

volScalarField dpore
(
    IOobject
    (
        "dpore",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar(porousProperties.lookup("dpore"), dimensionSet(0,1,0,0,0,0,0), porousProperties.lookup("dpore"))
);

volScalarField dpa
(
    IOobject
    (
        "dpa",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar(porousProperties.lookup("dpa"), dimensionSet(0,1,0,0,0,0,0), porousProperties.lookup("dpa"))
);

volScalarField dc
(
    IOobject
    (
        "dc",
        runTime.timeName(),
        mesh
   ),
    mesh,
   dimensionedScalar(porousProperties.lookup("dc"), dimensionSet(0,1,0,0,0,0,0), porousProperties.lookup("dc"))
);
    
volScalarField rhoC
(
    IOobject
    (
        "rhoC",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar(porousProperties.lookup("rhoC"), dimensionSet(1,-3,0,0,0,0,0),  porousProperties.lookup("rhoC"))
);

volScalarField cpC
(
    IOobject
    (
        "cpC",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar(porousProperties.lookup("cpC"), dimensionSet(0,2,-2,-1,0,0,0),  porousProperties.lookup("cpC"))
);

volScalarField kC
(
    IOobject
    (
        "kC",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar(porousProperties.lookup("kC"), dimensionSet(1,1,-3,-1,0,0,0),  porousProperties.lookup("kC"))
);

volScalarField e2
(
    IOobject
    (
        "e2",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("e2", dimensionSet(0,0,0,0,0,0,0), porosity2)
);

   
word bedPorosity(porousProperties.lookup("bedPorosity"));
    if ((bedPorosity != "experimental") && (bedPorosity != "calculated"))
    {
        FatalError << "in porousProperties, unknown bedPorosity type " << nl
            << "Valid types are: experimental or calculated."
            << bedPorosity << abort(FatalError);
    }

scalar porosity1 = 0.0;

   if(bedPorosity == "experimental")
   {

	porosity1 = readScalar(porousProperties.lookup("porosityb"));
    
   }
   else  //calculated
   {

	porosity1 = 0.368*(1+(2*(dpa[0]/dc[0])*((0.476/0.368)-1)));

    }

volScalarField e1
(
    IOobject
   (
        "e1",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("e1", dimensionSet(0,0,0,0,0,0,0), porosity1)
);


volScalarField et
(
    IOobject
    (
        "et",
        runTime.timeName(),
        mesh
    ),
    e1 + e2
);

volScalarField ei = e1/et;
volScalarField es = e2/et;

volScalarField as
(
    IOobject
    (
        "as",
        runTime.timeName(),
        mesh
    ),
    6*(1- e1)/(dpa*phi0)
);

Info << "eb=" << e1[0] << endl;
Info << "as=" << as[0] << endl;

volScalarField onebye1
(
    IOobject
    (
        "onebye1",
        runTime.timeName(),
        mesh
    ),
    1/e1
);

volScalarField onemenuse1
(
    IOobject
    (
        "onemenuse1",
        runTime.timeName(),
        mesh
    ),
    1 - e1
);

volScalarField onebye13
(
    IOobject
    (
        "onebye13",
        runTime.timeName(),
        mesh
    ),
    1/Foam::pow(e1,3)
);

volScalarField onemenuset
(
    IOobject
    (
        "onemenuset",
        runTime.timeName(),
        mesh
    ),
    1 - e1 - e2
);

volScalarField onemenuse2
(
    IOobject
    (
        "onemenuse2",
        runTime.timeName(),
        mesh
    ),
    1 - e2
);

volScalarField tau1
(
    IOobject
    (
        "tau1",
        runTime.timeName(),
        mesh
    ),
    (1.23*Foam::pow(1-e1, 1.33333))/(e1*Foam::pow(phi0, 2))
);

Info << "Bed tortuosity=" << tau1[0] << endl;

volScalarField tau2
(
    IOobject
    (
        "tau2",
        runTime.timeName(),
        mesh
    ),
    (1.23*Foam::pow(1-e2, 1.33333))/(e2*Foam::pow(phi0, 2))
);

Info << "Particle tortuosity=" << tau2[0] << endl;

volScalarField taut
(
    IOobject
    (
        "taut",
        runTime.timeName(),
        mesh
    ),
    tau1 + tau2
);

volScalarField EbyTau1
(
    IOobject
    (
        "EbyTau1",
        runTime.timeName(),
        mesh
    ),
    e1/tau1
);

volScalarField EbyTau2
(
    IOobject
    (
        "EbyTau2",
        runTime.timeName(),
        mesh
    ),
    e2/tau2
);

//Permeability of porous media

volScalarField Kp
(
   IOobject
    (
       "Kp",
       runTime.timeName(),
       mesh
   ),
   Foam::pow(dpa, 2)*Foam::pow(e1, 3)/(180*Foam::pow(onemenuse1, 2))
);

//Solving Darcy-Forschheimer terms in porous media

volScalarField Da
(
   IOobject
    (
       "Da",
       runTime.timeName(),
       mesh
   ),
   1/Kp
);

volScalarField e15 = Foam::pow(e1, 5);

volScalarField Fo
(
   IOobject
    (
       "Fo",
       runTime.timeName(),
       mesh
   ),
   1.8*e1/Foam::pow(180*e15*Kp, 0.5)
);

Info << "Permeability=" << Kp[0] << endl;
Info << "Da=" << Da[0] << endl;
Info << "Fo=" << Fo[0] << endl;
Info<< "\n" << endl;

//  necessaries units

volScalarField dimT
(
    IOobject
    (
        "dimT",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("dimT", dimensionSet(0,0,0,1,0,0,0), scalar(1))
);

volScalarField dimP
(
    IOobject
    (
        "dimP",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("dimP", dimensionSet(1,-1,-2,0,0,0,0), scalar(1))
);

volScalarField dimh
(
    IOobject
    (
        "dimh",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("dimh", dimensionSet(0,2,-2,0,0,0,0), scalar(1))
);

