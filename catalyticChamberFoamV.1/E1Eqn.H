{
    volScalarField& he1 = thermo1.he();
    volScalarField& hT1 = heatTransfer1;
    //volScalarField& hTg1 = heatTransferg1;
    //volScalarField& hTw = heatTransferw;
    //volScalarField& HmassT1 = sumHMassTransfer1;
    //volScalarField& massT = sumMassTransfer1;
    volScalarField& Mhe1 = Mdifhe1; //ok

    fvScalarMatrix E1Eqn
    (
        fvm::ddt(e1, rho1, he1)
      + mvConvection->fvmDiv(phi1, he1)
      + fvc::ddt(e1, rho1, K) 
      + fvc::div(phi1, K)
      + (
            he1.name() == "e"
          ? fvc::div
            (
                fvc::absolute(phi1/fvc::interpolate(rho1), U1),
                p1,
                "div(phiv,p)"
            )
          : -dpdt
        )
      - fvm::laplacian
        (
            fvc::interpolate(e1)
           *fvc::interpolate(turbulence->alphaEff()),
            he1
        )
     ==
        e1*Qdot1
      - hT1
      //- hTg1
      //- hTw
      //- fvm::Sp(massT, he1)
      //- HmassT1
      - Mhe1		//OK
      + fvOptions(rho1, he1)
    );

    E1Eqn.relax();

    fvOptions.constrain(E1Eqn);

    E1Eqn.solve();

    fvOptions.correct(he1);

    thermo1.correct();

    Info<< "min/max(T1) = "
        << min(T1).value() << ", " << max(T1).value() << endl;
}
