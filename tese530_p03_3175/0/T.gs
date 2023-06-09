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
    class       volScalarField;
    location    "0";
    object      T.stagnant;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 0 1 0 0 0];

internalField   uniform 530;

boundaryField
{
    inlet
    {
        type            fixedValue;
        value           uniform 530;
    }
   wall
    {
        type            fixedValue;
        value           uniform 530;
    }
    outlet
    {
        type             inletOutlet;
        phi		 phi2;
        inletValue       uniform 530;
        value            uniform 530;
    }
     eixo
    {
    type                 empty;
    }
    wedge1
    {
    type                wedge;
    }
    wedge2
    {
    type               wedge;
    }
}


// ************************************************************************* //
