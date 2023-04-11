IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties." + phase1Name,    // dictionary name
            runTime.constant(),     // dict is found in      "constant"
            mesh,                   // registry for the dict
            IOobject::MUST_READ,    // must exist, otherwise failure
            IOobject::NO_WRITE      // dict is only read by the solver
    )
    );
        
   PtrList<dimensionedScalar> sigma(Y1.size());
   PtrList<dimensionedScalar> epsilon(Y1.size());
   
   forAll(Y1, i)
    {
        word specieName = Y1[i].name();

        const dictionary& subDict2 = transportProperties.subDict(specieName);

        sigma.set
        (
             i,
             new dimensionedScalar(subDict2.lookup("sigma"))
        );
        
        epsilon.set
        (
             i,
             new dimensionedScalar(subDict2.lookup("epsilon"))
        );

    }
    
IOdictionary transportProperties2
    (
        IOobject
        (
            "transportProperties." + phase2Name,    // dictionary name
            runTime.constant(),     // dict is found in      "constant"
            mesh,                   // registry for the dict
            IOobject::MUST_READ,    // must exist, otherwise failure
            IOobject::NO_WRITE      // dict is only read by the solver
    )
    );
        
   PtrList<dimensionedScalar> sigma2(Y2.size());
   PtrList<dimensionedScalar> epsilon2(Y2.size());
   
   forAll(Y2, i)
    {
        word specieName2 = Y2[i].name();

        const dictionary& subDict2 = transportProperties2.subDict(specieName2);

        sigma2.set
        (
             i,
             new dimensionedScalar(subDict2.lookup("sigma"))
        );
        
        epsilon2.set
        (
             i,
             new dimensionedScalar(subDict2.lookup("epsilon"))
        );

    }
    