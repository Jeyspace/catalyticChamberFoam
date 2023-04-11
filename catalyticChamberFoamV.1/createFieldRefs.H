const volScalarField& psi1 = thermo1.psi();
const volScalarField& T1 = thermo1.T();
const volScalarField& T2 = thermo2.T();
const label inertIndex1(composition1.species()[inertSpecie1]);
const label inertIndex2(composition2.species()[inertSpecie2]);

scalar dtChem = refCast<const psiChemistryModel>(chemistry).deltaTChem()[0];

    PtrList<gasHThermoPhysics> specieData1(Y1.size());
    forAll(specieData1, i)
    {
        specieData1.set
        (
            i,
            new gasHThermoPhysics
            (
                dynamic_cast<const reactingMixture<gasHThermoPhysics>&>
                    (thermo1).speciesData()[i]
            )
        );
    }


    PtrList<gasHThermoPhysics> specieData2(Y2.size());
    forAll(specieData2, i)
    {
        specieData2.set
        (
            i,
            new gasHThermoPhysics
            (
                dynamic_cast<const reactingMixture<gasHThermoPhysics>&>
                    (thermo2).speciesData()[i]
            )
        );
    }


