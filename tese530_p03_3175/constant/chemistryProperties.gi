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
    object      chemistryProperties.interstitial;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

chemistryType
{
    chemistrySolver   EulerImplicit;
    chemistryThermo   psi;
}

chemistry           on;

initialChemicalTimeStep 1e-07;

EulerImplicitCoeffs
{
    cTauChem        1;
    equilibriumRateLimiter off;
}

odeCoeffs
{
    solver          Rosenbrock43;
    absTol          1e-12;
    relTol          0.01;
}

// ************************************************************************* //
