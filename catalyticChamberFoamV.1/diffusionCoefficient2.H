Info<< "Solving Mass Transfer 2 Model:\n" << endl;

Info<< "Solving SumWt2\n" << endl;

// Total sum of Mass: sumWt2=Y2[0]/mW0+Y2[1]/mW2+...+Y2[i]/mWi

volScalarField sumWt2 = 0.0*Y2[0];

forAll(Y2, i)
{   
    sumWt2 += Y2[i]/specieData2[i].W();
}
sumWt2.max(1e-2);

Info<< "sumWt2(max, min) =" << max(sumWt2).value() << ","<< min(sumWt2).value() <<endl;

volScalarField Wt2 = 1/sumWt2; //M=2/sumWt2

Info<< "Wt2(max, min) =" << max(Wt2).value() << ","<< min(Wt2).value() <<endl;

Info<< "Solving X2[i] and difX2[i]\n" << endl;

// definition of volScalarField 2 and molar fraction Xi=Yi*M/Mi

PtrList<volScalarField> X2(Y2.size()); 
PtrList<volScalarField> difX2(Y2.size()); 

forAll(Y2, i)
{ 

    volScalarField& Y2i = Y2[i];
    
     word nameX2("X2." + Y2[i].name());
    
            X2.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        nameX2,
                        runTime.timeName(),
                        mesh
                    ),
                    Y2i*Wt2/specieData2[i].W()
                )
            );

     word namedifX2("difX2." + Y2[i].name());
        
            difX2.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        namedifX2,
                        runTime.timeName(),
                        mesh
                    ),
                    1 - X2[i]
                )
            );

}


Info<< "Solving MDij2i and sumXDij2i\n" << endl;

Info<< "Solving fraction of T and p for the diffusion\n" << endl; //for difussion coefficient correction

volScalarField fT2
(
    IOobject
    (
        "fT2",
        runTime.timeName(),
        mesh
    ),
    (101325*dimP/p1)*Foam::pow(thermo2.T()/(273.15*dimT), 1.75)
);

//Mass diffusion coeficient calculation

PtrList<volScalarField> sumXDij2(Y2.size());  //sum for i-specie Xj*(2/Dij0)
 
forAll(Y2, i)
{ 
        
	volScalarField sumXDij20 = 0.0*Y2[0];
    
        forAll(Y2, j)
        {
            if(i!=j)
            {

                volScalarField& X2j = X2[j];
  
                sumXDij20 += X2j*MDij(i,j)/fT2;
                 
            }
            
        }
                   
            word nameSumXDij2("sumXDij2." + Y2[i].name());
        
            sumXDij2.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        nameSumXDij2,
                        runTime.timeName(),
                        mesh
                    ),
                    sumXDij20
                )
            );
        
}

Info<< "Solving Deff2\n" << endl;

volScalarField D2s
(
    IOobject
    (
        "D2s",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    EbyTau2*turbulence->muEff()/thermo2.rho()
);


//mixture effective diffusion
PtrList<volScalarField> Deff2(Y2.size());       //Efective diffusion diffusion of i-specie in mixture gas2

forAll(Y2, i)
{ 

	volScalarField L2 = 3.10609129e-4*T2/(Foam::pow(sigma[i].value(), 2)*p1);
	volScalarField Kn2 = L2/dpore;
	scalar Kn2max = max(Kn2).value();

	word nameDeff2("Deff2." + Y2[i].name());
                    
	Deff2.set
	(
		i,
		new volScalarField
		(
		IOobject
		(
		nameDeff2,
		runTime.timeName(),
		mesh,
        	IOobject::NO_READ,
        	IOobject::AUTO_WRITE
                             
		),
		D2s
                         
		)
                        
	);

	forAll(Deff2[i], cellI)
	{

            if(Y2[i][cellI] <= 0)
            {
		Deff2[i][cellI] = 1e-50;
            }
            else if(sumXDij2[i][cellI] <= 0)
            {
		Deff2[i][cellI] = 1e-50;
            }
            else if(Y2[i][cellI] == 1)
            {
		Deff2[i][cellI] = D2s[cellI];
            }
            else if(difX2[i][cellI] <= 0)
            {
		Deff2[i][cellI] = D2s[cellI];
            }
            else
            { 
                        
            if(Kn2[cellI] > 10)
            {
                        
            Deff2[i][cellI] = EbyTau2[cellI]*0.485*dpore[cellI]*Foam::pow(T2[cellI]/specieData2[i].W(), 0.5);

            }
            else if(Kn2[cellI] > 0.1)
            {

		Deff2[i][cellI] = 1/(1/(EbyTau2[cellI]*difX2[i][cellI]/sumXDij2[i][cellI])) + (1/(EbyTau2[cellI]*0.485*dpore[cellI]*Foam::pow(T2[cellI]/specieData2[i].W(), 0.5)));

            }
            else
            {
            
            Deff2[i][cellI] = EbyTau2[cellI]*difX2[i][cellI]/sumXDij2[i][cellI];

            }  
        
            }
	}

	if(Kn2max > 10)
	{
            
	Info << "Knudsen Diffusion for "  << Y2[i].name() << " : " << "Deff2(max, min): (" << max(Deff2[i]).value() << " , " << min(Deff2[i]).value() << ")" << endl;

	}
	else if(Kn2max > 0.1)
	{

	Info << "Knudsen + Mixture Diffusion for " << Y2[i].name() << " : " << "Deff2(max, min): (" << max(Deff2[i]).value() << " , " << min(Deff2[i]).value() << ")" << endl;
                       
	}
	else
	{
                   
	Info << "Mixture Diffusion for " << Y2[i].name() << " : " << "Deff2(max, min): (" << max(Deff2[i]).value() << " , " << min(Deff2[i]).value() << ")" << endl;

	}  

Deff2[i].max(0.0);
    
}

