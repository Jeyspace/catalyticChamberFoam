Info<< "Solving heatTransfer2\n" << endl;

Info<< "Solving T2-Tc\n" << endl;

volScalarField difT2
(
    IOobject
    (
        "difT2",
        runTime.timeName(),
        mesh
    ),
     thermo2.T() - Tc
);

Info << "T2-Tc (max) = " << max(difT2).value() << " ; "<< "T2-Tc (min) = " << min(difT2).value() << endl;

volScalarField heatTransfer2
(
    IOobject
    (
        "heatTransfer2",
        runTime.timeName(),
        mesh
    ),
    6*difT2*hh/dpa
);

//Info << "heatTransfer2(max) = " << max(heatTransfer2).value() << " , "<< "heatTransfer2(min) = " << min(heatTransfer2).value() << endl;

//Info<< "Solving heatTransferg2\n" << endl;

//volScalarField heatTransferg2
//(
//    IOobject
//    (
//        "heatTransferg2",
//        runTime.timeName(),
//        mesh
//    ),
//    3*difTg*hh/dpa
//);

//Info << "heatTransferg2(max) = " << max(heatTransferg2).value() << " ; "<< "heatTransferg2(min) = " << min(heatTransferg2).value() << endl;

//Info<< "Solving sumTmassTransfer Coefficient of mass difussion\n" << endl;

//scalar sumMassTransfer1 = 0.0;

//forAll(Y1, i)
//{    
//   sumMassTransfer1 += massTransfer1[i][0];
//}

//scalar sumMassTransfer2 = 0.0;

//forAll(Y2, i)
//{    
//    sumMassTransfer2 += massTransfer2[i][0];
//}

//scalar sumMassT = sumMassTransfer1 + sumMassTransfer2;


//volScalarField sumMassTT
//(
//    IOobject
//   (
//        "sumMassTT",
//        runTime.timeName(),
//       mesh
//    ),
//    mesh,
//    dimensionedScalar("sumMassTT",  dimensionSet(1,-3,-1,0,0,0,0), sumMassT)
//);


