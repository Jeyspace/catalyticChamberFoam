{

    volScalarField& hT3 = heatTransfer3;
    volScalarField& hT4 = heatTransfer4;

    fvScalarMatrix TcEqn
    (
        fvm::ddt(onemenuset, Tc)
      - fvm::laplacian
        (
	    fvc::interpolate(onemenuset)
           *fvc::interpolate(DTc),
            Tc
        )
     ==
	onemenuset*Qdot2/(rhoC*cpC)
      + hT3/(rhoC*cpC)
      + hT4/(rhoC*cpC)
      + fvOptions(Tc)
    );

    TcEqn.relax();

    fvOptions.constrain(TcEqn);

    TcEqn.solve();

    fvOptions.correct(Tc);

    Info<< "min/max(Tc) = "
        << min(Tc).value() << ", " << max(Tc).value() << endl;
}


