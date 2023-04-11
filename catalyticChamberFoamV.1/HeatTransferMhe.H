Info<< "Solving sumTmassTransfer Coefficient of mass difussion\n" << endl;

volScalarField sumMassTransfer1 = 0.0*massTransfer1[1];

forAll(Y1, i)
{
	sumMassTransfer1 += massTransfer1[i];
}

sumMassTransfer1.relax();

volScalarField sumMassTransfer2 = 0.0*massTransfer2[1];

forAll(Y2, i)
{     
	sumMassTransfer2 += massTransfer2[i];
}

sumMassTransfer2.relax();

volScalarField& he1 = thermo1.he();
volScalarField& he2 = thermo2.he();
volScalarField Mdifhe1 = sumMassTransfer1*(he1 - he2);
volScalarField Mdifhe2 = sumMassTransfer2*(he1 - he2);
