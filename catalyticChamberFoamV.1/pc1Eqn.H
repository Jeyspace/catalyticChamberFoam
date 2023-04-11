rhom = (e1*thermo1.rho() + e2*thermo2.rho())/et;

volScalarField rAU(1.0/U1Eqn.A());
volScalarField rAtU(1.0/(1.0/rAU - U1Eqn.H1()));
volVectorField HbyA(constrainHbyA(rAU*U1Eqn.H(), U1, p1));

if (pimple.nCorrPISO() <= 1)
{
    tU1Eqn.clear();
}

if (pimple.transonic())
{
    surfaceScalarField phid
    (
        "phid",
        fvc::interpolate(psi1)
       *(
            fvc::flux(HbyA)
          + fvc::interpolate(rhom*rAU)*fvc::ddtCorr(rhom, U1, phi1)
           /fvc::interpolate(rhom)
        )
    );


    MRF.makeRelative(fvc::interpolate(psi1), phid);

    surfaceScalarField phic
    (
        "phic",
        fvc::interpolate(rhom*(rAtU - rAU))*fvc::snGrad(p1)*mesh.magSf()
    );

    HbyA -= (rAU - rAtU)*fvc::grad(p1);

    volScalarField rhorAtU("rhorAtU", rhom*rAtU);

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p1Eqn
        (
            fvm::ddt(psi1, p1)
          + fvm::div(phid, p1)
          + fvc::div(phic)
          - fvm::laplacian(rhorAtU, p1)
         ==
            fvOptions(psi1, p1, rho1.name())
        );

        p1Eqn.solve(mesh.solver(p1.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phi1 == phic + p1Eqn.flux();
        }
    }
}
else
{
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        (
            fvc::flux(rhom*HbyA)
          + fvc::interpolate(rhom*rAU)*fvc::ddtCorr(rhom, U1, phi1)
        )
    );

    MRF.makeRelative(fvc::interpolate(rhom), phiHbyA);

    phiHbyA += fvc::interpolate(rhom*(rAtU - rAU))*fvc::snGrad(p1)*mesh.magSf();
    HbyA -= (rAU - rAtU)*fvc::grad(p1);

    volScalarField rhorAtU("rhorAtU", rhom*rAtU);

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p1, rho1, U1, phiHbyA, rhorAtU, MRF);

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p1Eqn
        (
            fvm::ddt(psi1, p1)
          + fvc::div(phiHbyA)
          - fvm::laplacian(rhorAtU, p1)
         ==
            fvOptions(psi1, p1, rho1.name())
        );

        p1Eqn.solve(mesh.solver(p1.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phi1 = phiHbyA + p1Eqn.flux();
        }
    }
}

#include "rho1Eqn.H"
#include "compressibleContinuityErrs.H"

// Explicitly relax pressure for momentum corrector
p1.relax();

U1 = HbyA - rAtU*fvc::grad(p1);
U1.correctBoundaryConditions();
fvOptions.correct(U1);
K = 0.5*magSqr(U1);

if (pressureControl.limit(p1))
{
    p1.correctBoundaryConditions();
    rhom = (e1*thermo1.rho() + e2*thermo2.rho())/et;
}

if (thermo1.dpdt())
{
    dpdt = fvc::ddt(p1);
}
