{
    volScalarField& he2 = thermo2.he();
    volScalarField& hT2 = heatTransfer2;
    //volScalarField& hTg2 = heatTransferg2;
    //volScalarField& HmassT2 = sumHMassTransfer2;
    //volScalarField& massT = sumMassTransfer2;
    volScalarField& Mhe2 = Mdifhe2; //ok

    fvScalarMatrix E2Eqn
    (
        fvm::ddt(e2, rho2, he2)
    ==
	e2*Qdot2
      - hT2
      //+ hTg2
      //+ fvm::Sp(massT, he2)
      //+ HmassT2
      + Mhe2		//OK
      + fvOptions(rho2, he2)
    );

   E2Eqn.relax();

   fvOptions.constrain(E2Eqn);

   E2Eqn.solve();

   fvOptions.correct(he2);

   thermo2.correct();

   Info<< "min/max(T2) = "
     << min(T2).value() << ", " << max(T2).value() << endl;
}


