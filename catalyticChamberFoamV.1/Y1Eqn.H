
tmp<fv::convectionScheme<scalar>> mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi1,
        mesh.divScheme("div(phi1,Y1i_h)")
    )
);

{
    reaction1->correct();
    Qdot1 = reaction1->Qdot();
    volScalarField Y1t(0.0*Y1[0]);

    forAll(Y1, i)
    {
        if (i != inertIndex1 && composition1.active(i))
        {
            volScalarField& Y1i = Y1[i];
            volScalarField& Deff1i = Deff1[i];
	    volScalarField& massTransfer1i = massTransfer1[i];
	    volScalarField Deffi = Defft1[i]; 

            fvScalarMatrix Y1iEqn
            (
               fvm::ddt(e1, rho1, Y1i)
             + mvConvection->fvmDiv(phi1, Y1i)
             - fvm::laplacian
               (
		fvc::interpolate(tau1)
               *fvc::interpolate(Deff1i)
               *fvc::interpolate(rho1), 
               Y1i
               ) 
             ==
               e1*reaction1->R(Y1i)
	     - massTransfer1i
             + fvOptions(rho1, Y1i)
            );

            Y1iEqn.relax();

            fvOptions.constrain(Y1iEqn);

            Y1iEqn.solve(mesh.solver("Y1i"));

            fvOptions.correct(Y1i);

            Y1i.max(0.0);

            Y1t += Y1i;
        }
    }

    Y1[inertIndex1] = scalar(1) - Y1t;
    Y1[inertIndex1].max(0.0);
}
