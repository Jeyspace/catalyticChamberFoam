
{
    Qdot2 = chemistry.Qdot()();
    volScalarField Y2t(0.0*Y2[0]);

    forAll(Y2, i)
    {
        if (i != inertIndex2 && composition2.active(i))
        {
            volScalarField& Y2i = Y2[i];
            volScalarField& massTransfer2i = massTransfer2[i];

            fvScalarMatrix Y2iEqn
            (
              fvm::ddt(e2, rho2, Y2i)	
             ==
	      chemistry.RR(i)
	    + massTransfer2i
            + fvOptions(rho2, Y2i)
            );
		
	    Y2iEqn.relax();

	    fvOptions.constrain(Y2iEqn);

            Y2iEqn.solve(mesh.solver("Y2i"));

            fvOptions.correct(Y2i);
         
            Y2i.max(0.0);

            Y2t += Y2i;
            
        }
    }

    Y2[inertIndex2] = /*Y1[inertIndex1]*/scalar(1) - Y2t;
    Y2[inertIndex2].max(0.0);

}
