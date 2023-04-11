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
    object      NH3.stagnant;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 0 0 0 0 0];

internalField   uniform 0;

boundaryField
{
   inlet
    {
        type            fixedValue;
        value           uniform 0;
    }
   wall
    {
        type            zeroGradient;
       
    }
    outlet
    {
        type            inletOutlet;
        phi             phi2;
        inletValue      $internalField;
        value           $internalField;
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
