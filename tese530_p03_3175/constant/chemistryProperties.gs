/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  4.x                                   |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant";
    object      chemistryProperties.stagnant;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

chemistryType
{
    chemistrySolver   /*ode*/EulerImplicit;
    chemistryThermo   psi;
}

chemistry       on;

initialChemicalTimeStep 1e-7;

EulerImplicitCoeffs
{
    cTauChem        1;
    equilibriumRateLimiter off;
}

odeCoeffs
{
    solver          /*seulex*/Rosenbrock43;
    absTol          1e-12;
    relTol          1e-1;
}


// ************************************************************************* //
