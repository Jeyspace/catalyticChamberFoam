// creating Phases throught words
    
    Info<< "Reading gasesNames.\n" << endl;
    IOdictionary phases
    (
        IOobject
        (
            "gases",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );
    
    const word phase1Name
    (
        phases.lookup("gas1")
    );
    
    const word phase2Name
    (
        phases.lookup("gas2")
    );





  

    

 
