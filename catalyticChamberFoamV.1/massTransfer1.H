Info<< "Solving Defft1 Coefficient of mass difussion\n" << endl;

PtrList<volScalarField> Defft1(Y1.size());   

forAll(Y1, i)
{
	volScalarField& Deff1i = Deff1[i];
	volScalarField& Deff2i = Deff2[i];

	volScalarField mDeffi = taut*Deff1i*Deff2i;
	volScalarField sDeffi = et*(Deff1i+Deff2i);

	if(mDeffi[0] <= 0)
	{
	
	word nameDefft1("Defft1." + Y1[i].name());
        
        Defft1.set
        (
            i,
         new volScalarField
         (
             IOobject
             (
                 nameDefft1,
              runTime.timeName(),
              mesh
                
            ),
         D1s*0
         )
         );

	}
	else if(sDeffi[0] <= 0)
	{

	word nameDefft1("Defft1." + Y1[i].name());
        
        Defft1.set
        (
            i,
         new volScalarField
         (
             IOobject
             (
                 nameDefft1,
              runTime.timeName(),
              mesh
                
            ),
         D1s*0
         )
         );

	}
	else
	{
	
	word nameDefft1("Defft1." + Y1[i].name());
        
        Defft1.set
        (
            i,
         new volScalarField
         (
             IOobject
             (
                 nameDefft1,
              runTime.timeName(),
              mesh
                
            ),
            mDeffi/sDeffi
         )
         );
	}

}

Info<< "Solving kc1 Coefficient of mass difussion\n" << endl;

PtrList<volScalarField> kc1(Y1.size());   

forAll(Y1, i)
{
	volScalarField Deffi = Defft1[i];        
        volScalarField Sc1 = mum/(Deffi*rhom);
        
        volScalarField Sh1 = scalar(2) + 1.1*Rep*cbrt(Sc1);
        
        word namekc1("kc1." + Y1[i].name());
        
        kc1.set
        (
            i,
         new volScalarField
         (
             IOobject
             (
                 namekc1,
              runTime.timeName(),
              mesh
                
            ),
             Deffi*Sh1/dpa
         )
            
        );    

	forAll(kc1[i], cellI)
	{
		if(Defft1[i][cellI] >= 1e-3)
		{
		kc1[i][cellI] = 0;
		}
	}

}

Info<< "Solving difY1\n" << endl;

PtrList<volScalarField> difY1(Y1.size()); 

forAll(Y1, i)
{

     volScalarField& Y1i = Y1[i];
     volScalarField& Y2i = Y2[i];
        
        word namedifY1("difY1." + Y1[i].name());
            
            difY1.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        namedifY1,
                        runTime.timeName(),
                        mesh
                    ),
                     rhom*(Y1i - Y2i)
                )
            );
}

Info<< "Solving massTransfer1\n" << endl;

PtrList<volScalarField> massTransfer1(Y1.size()); 

forAll(Y1, i)
{
        volScalarField& kc1i = kc1[i];
        volScalarField& difY1i = difY1[i];
        
        word nameMassTransfer1("massTransfer1." + Y1[i].name());
            
            massTransfer1.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        nameMassTransfer1,
                        runTime.timeName(),
                        mesh                     
                    ),
                    as*kc1i*difY1i
                )
            );
            

}

