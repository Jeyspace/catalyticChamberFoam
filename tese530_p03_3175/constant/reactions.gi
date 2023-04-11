/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  5.x                                   |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       IOobject;
    location    "constant";
    object      reactions.gi;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

elements
(
N
H
);

species
(
    N2H4
    H2
    NH3
    N2
    Air
);

reactions
{
    HomoHydrazineReaction
    {
        type     irreversibleArrheniusReaction;
        reaction "2N2H4 = 2NH3 + N2 + H2";
        A        2.4e10;
        beta     0;
        Ta       18333.33;
    }
}
