Info<< "Solving T1-Tc\n" << endl;

volScalarField difT1
(
    IOobject
    (
        "difT1",
        runTime.timeName(),
        mesh
    ),
    thermo1.T() - Tc
);


Info << "T1-Tc (max) = " << max(difT1).value() << " ; "<< "T1-Tc (min) = " << min(difT1).value() << endl;

Info<< "Solving heatTransfer1\n" << endl;

volScalarField heatTransfer1
(
    IOobject
    (
        "heatTransfer1",
        runTime.timeName(),
        mesh
    ),
    /*0.5*/as*difT1*hh
);

Info << "heatTransfer1(max) = " << max(heatTransfer1).value() << " ; "<< "heatTransfer1(min) = " << min(heatTransfer1).value() << endl;

//Info<< "Solving T1-T2\n" << endl;

//volScalarField difTg
//(
 //   IOobject
 //   (
 //       "difTg",
 //       runTime.timeName(),
  //      mesh
 //   ),
 //   thermo1.T() - thermo2.T()
//);


//Info << "T1-T2 (max) = " << max(difTg).value() << " ; "<< "T1-T2 (min) = " << min(difTg).value() << endl;

//Info<< "Solving heatTransferg1\n" << endl;

//volScalarField heatTransferg1
//(
//    IOobject
//    (
//        "heatTransferg1",
//        runTime.timeName(),
//        mesh
//    ),
//    0.5*as*difTg*hh
//);

//Info << "heatTransferg1(max) = " << max(heatTransferg1).value() << " ; "<< "heatTransferg1(min) = " << min(heatTransferg1).value() << endl;
