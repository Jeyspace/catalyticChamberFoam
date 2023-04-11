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
    class       volScalarField;
    location    "0";
    object      Ydefault.gs;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 0 0 0 0 0];

internalField   uniform 0;

boundaryField
{
    outlet
    {
        type            calculated;
        value           uniform 0;
    }
    inlet
    {
        type            calculated;
        value           uniform 0;
    }
    wall
    {
        type            calculated;
        value           uniform 0;
    }
    wedge1
    {
        type            wedge;
    }
    wedge2
    {
        type            wedge;
    }
    eixo
    {
        type            empty;
    }
}


// ************************************************************************* //
