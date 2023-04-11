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
    object      transportProperties.gi; //REF.A SOFTWARE PACKAGE FOR THE EVALUATION OF GAS-PHASE, MULTICOMPONENT TRANSPORT PROPERTIES
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


//Transport coefficient Type: Wakao or Bird(Kesten used)
//transport	Wakao;

N2H4.gi
{
sigma    sigma [0 0 0 1 0 0 0]4.230;
epsilon	 epsilon [0 0 0 0 0 0 0]205;
}

H2.gi
{
sigma    sigma [0 0 0 1 0 0 0]2.92;
epsilon	 epsilon [0 0 0 0 0 0 0]38.000;
}

NH3.gi
{
sigma    sigma [0 0 0 1 0 0 0]2.920;
epsilon	 epsilon [0 0 0 0 0 0 0]481;
}

N2.gi
{
sigma    sigma [0 0 0 1 0 0 0]3.621;
epsilon	 epsilon [0 0 0 0 0 0 0]97.530;
}

Air.gi
{
sigma    sigma [0 0 0 1 0 0 0]3.617;
epsilon	 epsilon [0 0 0 0 0 0 0]97;
}


