//Info<< "Solving Defft2  Coefficient of mass difussion 2\n" << endl;

Info<< "Solving massTransfer2\n" << endl;

PtrList<volScalarField> massTransfer2(Y2.size()); 

forAll(Y2, i)
{
        volScalarField& kc2i = kc1[i];
        volScalarField& difY2i = difY1[i];
        
        word nameMassTransfer2("massTransfer2." + Y2[i].name());
            
            massTransfer2.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        nameMassTransfer2,
                        runTime.timeName(),
                        mesh
                    ),
                    6*kc2i*difY2i/dpa
                )
            );

}





