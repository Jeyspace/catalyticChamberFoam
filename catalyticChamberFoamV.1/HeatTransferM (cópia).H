Info<< "Solving species enthalpy Coefficient of mass difussion\n" << endl;

PtrList<volScalarField> hei1(Y1.size());  

forAll(Y1, i)
{

	scalar T1i = 273.15;
	scalar p1i = 101325;

	scalar hei10 = specieData1[i].he(p1i, T1i);

        word namehei1("hei1." + Y1[i].name());
        
        hei1.set
        (
            i,
         new volScalarField
         (
             IOobject
             (
                 namehei1,
              runTime.timeName(),
              mesh
                
            ),
	        dimh*hei10
         )
         );

}

forAll(Y1, i)
{
	volScalarField& hei1i = hei1[i];

	forAll(hei1i, cellI)
	{

	scalar T1i = T1[cellI];
	scalar p1i = p1[cellI];

	hei1i[cellI] = specieData1[i].he(p1i, T1i);

	}

}


PtrList<volScalarField> hei2(Y1.size());  

forAll(Y2, i)
{

	scalar T2i = 273.15;
	scalar p1i = 101325;

	scalar hei20 = specieData2[i].he(p1i, T2i);

        word namehei2("hei2." + Y2[i].name());
        
        hei2.set
        (
            i,
         new volScalarField
         (
             IOobject
             (
                 namehei2,
              runTime.timeName(),
              mesh
                
            ),
	        dimh*hei20
         )
         );

}

forAll(Y2, i)
{
	volScalarField& hei2i = hei2[i];

	forAll(hei2i, cellI)
	{

	scalar T2i = T2[cellI];


	scalar p1i = p1[cellI];

	hei2i[cellI] = specieData2[i].he(p1i, T2i);

	}

}

Info<< "Solving difhi\n" << endl;

PtrList<volScalarField> difhi(Y1.size()); 

forAll(Y1, i)
{
     volScalarField& hei1i = hei1[i];
     volScalarField& hei2i = hei2[i];
        
        word namedifhi("difhi." + Y1[i].name());
            
            difhi.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        namedifhi,
                        runTime.timeName(),
                        mesh
                    ),
                    thermo1.rho()*hei1i - thermo2.rho()*hei2i
                )
            );
            
}

Info<< "Solving HmassTransfer1\n" << endl;

PtrList<volScalarField> HmassTransfer1(Y1.size()); 

forAll(Y1, i)
{

        volScalarField& kc1i = kc1[i];
        volScalarField& difhii = difhi[i];
        
        word nameHMassTransfer1("HmassTransfer1." + Y1[i].name());
            
            HmassTransfer1.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        nameHMassTransfer1,
                        runTime.timeName(),
                        mesh                     
                    ),
                    as*kc1i*difhii
                )
            );

}

Info<< "Solving HmassTransfer2\n" << endl;

PtrList<volScalarField> HmassTransfer2(Y2.size()); 

forAll(Y2, i)
{

        volScalarField& kc1i = kc1[i];
        volScalarField& difhii = difhi[i];
        
        word nameHMassTransfer2("HmassTransfer2." + Y1[i].name());
            
            HmassTransfer2.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        nameHMassTransfer2,
                        runTime.timeName(),
                        mesh                     
                    ),
                    6*kc1i*difhii/dpa
		)
            );

}

Info<< "Solving sumHmassTransfer1 Coefficient of mass difussion\n" << endl;

volScalarField sumHMassTransfer1 = 0.0*HmassTransfer1[0];

forAll(Y1, i)
{
	sumHMassTransfer1 += HmassTransfer1[i];
}

//sumHMassTransfer1.relax();

Info<< "Solving sumHmassTransfer2 Coefficient of mass difussion\n" << endl;

volScalarField sumHMassTransfer2 = 0.0*HmassTransfer2[0];

forAll(Y2, i)
{    
	sumHMassTransfer2 += HmassTransfer2[i];
}

//sumHMassTransfer2.relax();

Info<< "Solving sumTmassTransfer Coefficient of mass difussion\n" << endl;

volScalarField sumMassTransfer1 = 0.0*massTransfer1[1];

forAll(Y1, i)
{
	sumMassTransfer1 += massTransfer1[i];
}

//sumMassTransfer1.relax();

volScalarField sumMassTransfer2 = 0.0*massTransfer2[1];

forAll(Y2, i)
{     
	sumMassTransfer2 += massTransfer2[i];
}

//sumMassTransfer2.relax();

	volScalarField& he1 = thermo1.he();
	volScalarField& he2 = thermo2.he();
	//Q = 2*mag(U1)*thermo1.rho()/dc;
	volScalarField Mdifhe = sumMassTransfer1*(he1 - he2);
	//volScalarField Qdifhe = Q*(he1 - he2);
