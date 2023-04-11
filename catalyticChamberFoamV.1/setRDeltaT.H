/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

{
    volScalarField& rDeltaT = trDeltaT.ref();

    const dictionary& pimpleDict = pimple.dict();

    // Maximum flow Courant number
    scalar maxCo(readScalar(pimpleDict.lookup("maxCo")));

    // Maximum time scale
    scalar maxDeltaT(pimpleDict.lookupOrDefault<scalar>("maxDeltaT", GREAT));

    // Smoothing parameter (0-1) when smoothing iterations > 0
    scalar rDeltaTSmoothingCoeff
    (
        pimpleDict.lookupOrDefault<scalar>("rDeltaTSmoothingCoeff", 0.1)
    );

    // Damping coefficient (1-0)
    scalar rDeltaTDampingCoeff
    (
        pimpleDict.lookupOrDefault<scalar>("rDeltaTDampingCoeff", 1.0)
    );

    // Maximum change in cell temperature per iteration
    // (relative to previous value)
    scalar alphaTemp(pimpleDict.lookupOrDefault("alphaTemp", 0.05));

    // Maximum change in cell concentration per iteration
    // (relative to reference value)
    scalar alphaY1(pimpleDict.lookupOrDefault("alphaY1", 1.0));

    Info<< "Time scales min/max:" << endl;

    // Cache old reciprocal time scale field
    volScalarField rDeltaT0("rDeltaT0", rDeltaT);

    // Flow time scale
    {
        rDeltaT.ref() =
        (
            fvc::surfaceSum(mag(phi1))()()
           /((2*maxCo)*mesh.V()*rho1())
        );

        // Limit the largest time scale
        rDeltaT.max(1/maxDeltaT);

        Info<< "    Flow        = "
            << 1/gMax(rDeltaT.primitiveField()) << ", "
            << 1/gMin(rDeltaT.primitiveField()) << endl;
    }

    // Heat release rate time scale
    if (alphaTemp < 1)
    {
        volScalarField::Internal rDeltaTT
        (
            mag(Qdot1)/(alphaTemp*rho1*thermo1.Cp()*T1)
        );

        Info<< "    Temperature = "
            << 1/(gMax(rDeltaTT.field()) + VSMALL) << ", "
            << 1/(gMin(rDeltaTT.field()) + VSMALL) << endl;

        rDeltaT.ref() = max(rDeltaT(), rDeltaTT);
    }

    // Reaction rate time scale
    if (alphaY1 < 1)
    {
        dictionary Y1ref(pimpleDict.subDict("Y1ref"));

        volScalarField::Internal rDeltaTY1
        (
            IOobject
            (
                "rDeltaTY1",
                runTime.timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar("rDeltaTY1", rDeltaT.dimensions(), 0)
        );

        bool foundY1 = false;

        forAll(Y1, i)
        {
            if (i != inertIndex1 && composition1.active(i))
            {
                volScalarField& Y1i = Y1[i];

                if (Y1ref.found(Y1i.name()))
                {
                    foundY1 = true;
                    scalar Y1refi = readScalar(Y1ref.lookup(Y1i.name()));

                    rDeltaTY1.field() = max
                    (
                        mag
                        (
                            reaction1->R(Y1i)().source()
                           /((Y1refi*alphaY1)*(rho1*mesh.V()))
                        ),
                        rDeltaTY1
                    );
                }
            }
        }

        if (foundY1)
        {
            Info<< "    Composition = "
                << 1/(gMax(rDeltaTY1.field()) + VSMALL) << ", "
                << 1/(gMin(rDeltaTY1.field()) + VSMALL) << endl;

            rDeltaT.ref() = max(rDeltaT(), rDeltaTY1);
        }
        else
        {
            IOWarningIn(args.executable().c_str(), Y1ref)
                << "Cannot find any active species in Y1ref " << Y1ref
                << endl;
        }
    }

    // Update tho boundary values of the reciprocal time-step
    rDeltaT.correctBoundaryConditions();

    // Spatially smooth the time scale field
    if (rDeltaTSmoothingCoeff < 1)
    {
        fvc::smooth(rDeltaT, rDeltaTSmoothingCoeff);
    }

    // Limit rate of change of time scale
    // - reduce as much as required
    // - only increase at a fraction of old time scale
    if
    (
        rDeltaTDampingCoeff < 1
     && runTime.timeIndex() > runTime.startTimeIndex() + 1
    )
    {
        rDeltaT = max
        (
            rDeltaT,
            (scalar(1) - rDeltaTDampingCoeff)*rDeltaT0
        );
    }

    // Update tho boundary values of the reciprocal time-step
    rDeltaT.correctBoundaryConditions();

    Info<< "    Overall     = "
        << 1/gMax(rDeltaT.primitiveField())
        << ", " << 1/gMin(rDeltaT.primitiveField()) << endl;
}


// ************************************************************************* //
