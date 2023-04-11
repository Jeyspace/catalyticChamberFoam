Info<< "Solving Mass Transfer 1 Model:\n" << endl;

Info<< "Solving SumWt1\n" << endl;

// Total sum of Mass: sumWt1=Y1[0]/mW0+Y1[1]/mW1+...+Y1[i]/mWi

volScalarField sumWt1 = 0.0*Y1[0];

forAll(Y1, i)
{   
    sumWt1 += Y1[i]/specieData1[i].W();
}

sumWt1.max(1e-2);

Info<< "sumWt1(max, min) =" << max(sumWt1).value() << ","<< min(sumWt1).value() <<endl;

volScalarField Wt1 = 1/sumWt1; //M=1/sumWt1

Info<< "Wt1(max, min) =" << max(Wt1).value() << ","<< min(Wt1).value() <<endl;

Info<< "Solving X1[i] and difX1[i]\n" << endl;

// definition of volScalarField 1 and molar fraction Xi=Yi*M/Mi

PtrList<volScalarField> X1(Y1.size()); 
PtrList<volScalarField> difX1(Y1.size()); 

forAll(Y1, i)
{ 

    volScalarField& Y1i = Y1[i];
    
     word nameX1("X1." + Y1[i].name());
    
            X1.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        nameX1,
                        runTime.timeName(),
                        mesh
                    ),
                    Y1i*Wt1/specieData1[i].W()
                )
            );

     word namedifX1("difX1." + Y1[i].name());
        
            difX1.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        namedifX1,
                        runTime.timeName(),
                        mesh
                    ),
                    1 - X1[i]
                )
            );

}


Info<< "Solving MDij1i and sumXDij1i\n" << endl;

//Mass diffusion coeficient calculation

Info<< "Solving fraction of T and p for the diffusion\n" << endl; //for difussion coefficient correction

volScalarField fT1
(
    IOobject
    (
        "fT1",
        runTime.timeName(),
        mesh
    ),
    (101325*dimP/p1)*Foam::pow(thermo1.T()/(273.15*dimT), 1.75)
);

PtrList<volScalarField> sumXDij1(Y1.size());  //sum for i-specie Xj*(1/Dij0)
 
forAll(Y1, i)
{ 
        
	volScalarField sumXDij10 = 0.0*Y1[0];
    
        forAll(Y1, j)
        {
            if(i!=j)
            {

                volScalarField& X1j = X1[j];
  
                sumXDij10 += X1j*MDij(i,j)/fT1;
                 
            }
            
        }
                   
            word nameSumXDij1("sumXDij1." + Y1[i].name());
        
            sumXDij1.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        nameSumXDij1,
                        runTime.timeName(),
                        mesh
                    ),
                    sumXDij10
                )
            );
        
}

Info<< "Solving Deff1\n" << endl;

volScalarField D1s
(
    IOobject
    (
        "D1s",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    EbyTau1*turbulence->muEff()/thermo1.rho()
);

volScalarField& T1 = thermo1.T();
volScalarField& p1 = thermo1.p();

//mixture effective diffusion
PtrList<volScalarField> Deff1(Y1.size());       //Efective diffusion diffusion of i-specie in mixture gas1

forAll(Y1, i)
{ 
	volScalarField L1 = 3.10609129e-4*T1/(Foam::pow(sigma[i].value(), 2)*p1);
	volScalarField Kn1 = L1/dpore;
	scalar Kn1max = max(Kn1).value();

	word nameDeff1("Deff1." + Y1[i].name());
                    
	Deff1.set
	(
		i,
		new volScalarField
		(
		IOobject
		(
		nameDeff1,
		runTime.timeName(),
		mesh,
        	IOobject::NO_READ,
        	IOobject::AUTO_WRITE
                             
		),
		D1s
                         
		)
                        
	);

	forAll(Deff1[i], cellI)
	{
            if(Y1[i][cellI] <= 0)
            {
		Deff1[i][cellI] = 1e-50;
            }
            else if(sumXDij1[i][cellI] <= 0)
            {
		Deff1[i][cellI] = 1e-50;
            }
            else if(Y1[i][cellI] == 1)
            {
		Deff1[i][cellI] = D1s[cellI];
            }
            else if(difX1[i][cellI] <= 0)
            {
		Deff1[i][cellI] = D1s[cellI];
            }
            else
            { 
                        
            if(Kn1[cellI] > 10)
            {
                        
            Deff1[i][cellI] = EbyTau1[cellI]*0.485*dpore[cellI]*Foam::pow(T1[cellI]/specieData1[i].W(), 0.5);

            }
            else if(Kn1[cellI] > 0.1)
            {

		Deff1[i][cellI] = 1/(1/(EbyTau1[cellI]*difX1[i][cellI]/sumXDij1[i][cellI])) + (1/(EbyTau1[cellI]*0.485*dpore[cellI]*Foam::pow(T1[cellI]/specieData1[i].W(), 0.5)));

            }
            else
            {
            
            Deff1[i][cellI] = EbyTau1[cellI]*difX1[i][cellI]/sumXDij1[i][cellI];

            }  
        
            }
	}

	if(Kn1max > 10)
	{
            
	Info << "Knudsen Diffusion for "  << Y1[i].name() << " : " << "Deff1(max, min): (" << max(Deff1[i]).value() << " , " << min(Deff1[i]).value() << ")" << endl;

	}
	else if(Kn1max > 0.1)
	{

	Info << "Knudsen + Mixture Diffusion for " << Y1[i].name() << " : " << "Deff1(max, min): (" << max(Deff1[i]).value() << " , " << min(Deff1[i]).value() << ")" << endl;
                       
	}
	else
	{
                   
	Info << "Mixture Diffusion for " << Y1[i].name() << " : " << "Deff1(max, min): (" << max(Deff1[i]).value() << " , " << min(Deff1[i]).value() << ")" << endl;

	}  

Deff1[i].max(0.0);
    
}

