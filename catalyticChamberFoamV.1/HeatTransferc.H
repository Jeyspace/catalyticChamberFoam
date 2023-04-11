
Info<< "Solving heatTransfer3\n" << endl;

volScalarField heatTransfer3
(
    IOobject
    (
        "heatTransfer3",
        runTime.timeName(),
        mesh
    ),
    6*difT1*hh/dpa
);

volScalarField heatTransfer4
(
    IOobject
    (
        "heatTransfer4",
        runTime.timeName(),
        mesh
    ),
    6*difT2*hh/dpa
);

