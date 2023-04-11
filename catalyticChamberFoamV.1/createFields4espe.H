#include "createRDeltaT.H"

#include "createPhases.H"	//Declaration of the phases

//Thermodinamyc model for Phase 1

Info<< "Creating reaction model for phase1Name\n" << endl;

autoPtr<combustionModels::psiCombustionModel> reaction1
(
    combustionModels::psiCombustionModel::New(mesh, "combustionProperties", phase1Name)
);

psiReactionThermo& thermo1 = reaction1->thermo();
thermo1.validate(args.executable(), "h", "e");

basicMultiComponentMixture& composition1 = thermo1.composition();
PtrList<volScalarField>& Y1 = composition1.Y();

const word inertSpecie1(thermo1.lookup("inertSpecie"));
if (!composition1.species().found(inertSpecie1))
{
    FatalIOErrorIn(args.executable().c_str(), thermo1)
        << "Inert specie " << inertSpecie1 << " not found in available species "
        << composition1.species()
        << exit(FatalIOError);
}

//Thermodinamyc model for Phase 2

// write base thermo fields - not registered since will be re-read by
// thermo package

Info<< "Creating base fields for time " << runTime.timeName() << endl;

{
  
    volScalarField Ydefault2
    (
        IOobject
        (
            IOobject::groupName("Ydefault", phase2Name),
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        ),
        mesh,
        dimensionedScalar("Ydefault", dimless, 0)
    );

    Ydefault2.write();
}

Info<< nl << "Reading thermophysicalProperties for Phase2Name" << endl;

autoPtr<psiChemistryModel> pChemistry(psiChemistryModel::New(mesh, phase2Name));

psiChemistryModel& chemistry = pChemistry();

psiReactionThermo& thermo2 = chemistry.thermo();
thermo2.validate(args.executable(), "h");

basicMultiComponentMixture& composition2 = thermo2.composition();
PtrList<volScalarField>& Y2 = composition2.Y();

const word inertSpecie2(thermo2.lookup("inertSpecie"));
if (!composition2.species().found(inertSpecie2))
{
   FatalIOErrorIn(args.executable().c_str(), thermo1)
        << "Inert specie " << inertSpecie2 << " not found in available species "
        << composition1.species()
        << exit(FatalIOError);
}

Info<< "Reading field Tc\n" << endl;

  volScalarField Tc
  (
       IOobject
      (
          "Tc",
          runTime.timeName(),
           mesh,
          IOobject::NO_READ,
          IOobject::AUTO_WRITE
      ),
      thermo2.T()
   );

Info<< "Creating Porous Properties\n" << endl;

#include "createPorousProperties.H"

// Density of gas mixture 1 and gas mixture 2 

volScalarField rho1
(
    IOobject
    (
        "rho1",
         runTime.timeName(),
         mesh
    ),
    thermo1.rho()
);

volScalarField rho2
(
    IOobject
    (
        "rho2",
         runTime.timeName(),
         mesh
    ),
    thermo2.rho()
);

volScalarField rhom
(
    IOobject
    (
        "rhom",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    (e1*thermo1.rho() + e2*thermo2.rho())/et
);

// velocity of phase U1

Info<< "Reading field U1\n" << endl;

volVectorField U1
(
    IOobject
    (
        IOobject::groupName("U", phase1Name),
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh 
);

volScalarField G
(
    IOobject
    (
        "G",
         runTime.timeName(),
         mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mag(U1)*thermo1.rho()
);

// velocity of phase2 U2 = 0, stagnant film
    
volVectorField U2
(
    IOobject
    (
        IOobject::groupName("U", phase2Name),
        runTime.timeName(),
        runTime,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("zero", dimVelocity, Zero)
);

volScalarField& p1 = thermo1.p();

#include "compressibleCreatePhi.H"	//phase1
#include "createPhi.H"			//phase2

// dates for output grafics

OFstream post(args.path()/"chemFoam.out1");
 post<< "# Time" << token::TAB 
        << "Temperature.catalyst [K]" << token::TAB
        << "Y1.catalyst" << token::TAB
        << "Y2.catalyst" << token::TAB
        << "Y3.catalyst" << token::TAB
        << "Y4.catalyst" << token::TAB
//        << "Y5.catalyst" << token::TAB 
        << "Temperature.gas [K]" << token::TAB
        << "Pressure.gas [Pa]" << token::TAB
        << "Y1.gas" << token::TAB
        << "Y2.gas" << token::TAB
        << "Y3.gas" << token::TAB
        << "Y4.gas" << token::TAB
//        << "Y5.gas" << token::TAB 
        << "Tsolid" << token::TAB
        << "Deff.catalyst" << token::TAB
        << "Deff.catalyst" << token::TAB
        << "Deff.catalyst" << token::TAB
        << "Deff.catalyst" << token::TAB
 //       << "Deff.catalyst" << token::TAB 
        << "Deff.gas" << token::TAB
        << "Deff.gas" << token::TAB
        << "Deff.gas" << token::TAB
        << "Deff.gas" << token::TAB
//        << "Deff.gas" << token::TAB 
        << "D2s" << token::TAB
	<< "D1s" << token::TAB 
        << "massTransfer11.catalyst" << token::TAB
        << "massTransfer12.catalyst" << token::TAB
        << "massTransfer13.catalyst" << token::TAB
        << "massTransfer14.catalyst" << token::TAB
  //      << "massTransfer15.catalyst" << token::TAB 
        << "massTransfer21.catalyst" << token::TAB
        << "massTransfer22.catalyst" << token::TAB
        << "massTransfer23.catalyst" << token::TAB
        << "massTransfer24.catalyst" << token::TAB << token::TAB << endl;
 //       << "massTransfer25.catalyst" << token::TAB  << token::TAB << endl;

 OFstream post2(args.path()/"chemFoam.out2");
 post2<< "# Time" << token::TAB 
        << "Temperature.catalyst [K]" << token::TAB
        << "Y1.catalyst" << token::TAB
        << "Y2.catalyst" << token::TAB
        << "Y3.catalyst" << token::TAB
        << "Y4.catalyst" << token::TAB
//        << "Y5.catalyst" << token::TAB 
        << "Temperature.gas [K]" << token::TAB
        << "Pressure.gas [Pa]" << token::TAB
        << "Y1.gas" << token::TAB
        << "Y2.gas" << token::TAB
        << "Y3.gas" << token::TAB
        << "Y4.gas" << token::TAB
//        << "Y5.gas" << token::TAB 
        << "Tsolid" << token::TAB
        << "Deff.catalyst" << token::TAB
        << "Deff.catalyst" << token::TAB
        << "Deff.catalyst" << token::TAB
        << "Deff.catalyst" << token::TAB
//        << "Deff.catalyst" << token::TAB 
        << "Deff.gas" << token::TAB
        << "Deff.gas" << token::TAB
        << "Deff.gas" << token::TAB
        << "Deff.gas" << token::TAB
//        << "Deff.gas" << token::TAB 
        << "D2s" << token::TAB
	<< "D2s" << token::TAB 
        << "massTransfer11.catalyst" << token::TAB
        << "massTransfer12.catalyst" << token::TAB
        << "massTransfer13.catalyst" << token::TAB
        << "massTransfer14.catalyst" << token::TAB
//        << "massTransfer15.catalyst" << token::TAB 
        << "massTransfer21.catalyst" << token::TAB
        << "massTransfer22.catalyst" << token::TAB
        << "massTransfer23.catalyst" << token::TAB
        << "massTransfer24.catalyst" << token::TAB << token::TAB << endl;
//        << "massTransfer25.catalyst" << token::TAB  << token::TAB << endl;

pressureControl pressureControl(p1, rho1, pimple.dict(), false);

mesh.setFluxRequired(p1.name());

Info << "Creating turbulence model.\n" << nl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho1,
        U1,
        phi1,
        thermo1
    )
);

// Set the turbulence into the reaction model
reaction1->setTurbulence(turbulence());


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
    dimensionedScalar("dpdt", p1.dimensions()/dimTime, 0)
);

Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K("K", 0.5*magSqr(U1));


multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

forAll(Y1, i)
{
    fields.add(Y1[i]);
}
fields.add(thermo1.he());

forAll(Y2, i)
{
    fields.add(Y2[i]);
}
fields.add(thermo2.he());

volScalarField Qdot1
(
    IOobject
    (
        "Qdot1",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Qdot1", dimEnergy/dimVolume/dimTime, 0.0)
);

volScalarField Qdot2
(
    IOobject
    (
        "Qdot2",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Qdot2", dimEnergy/dimVolume/dimTime, 0.0)
);

#include "createMRF.H"


