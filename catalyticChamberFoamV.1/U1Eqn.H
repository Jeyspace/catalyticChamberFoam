// Solve the Momentum equation

MRF.correctBoundaryVelocity(U1);
//volScalarField& massT = sumMassTransfer1;

tmp<fvVectorMatrix> tU1Eqn
(
    fvm::ddt(onebye1, rho1, U1)   
  + onebye13*fvm::div(phi1, U1)
  + onebye1*MRF.DDt(rho1, U1)
  + turbulence->divDevRhoReff(U1)
  - fvm::laplacian
    (
      fvc::interpolate(onebye1)
     *fvc::interpolate(thermo1.mu()),
      U1
    )
 ==
  - fvm::Sp(turbulence->muEff()*Da, U1)
  - fvm::Sp(Fo*thermo1.rho()*mag(U1), U1)
  //- fvm::Sp(massT, U1)
  + fvOptions(rho1, U1)
);

fvVectorMatrix& U1Eqn = tU1Eqn.ref();

U1Eqn.relax();

fvOptions.constrain(U1Eqn);

if (pimple.momentumPredictor())
{
    solve(U1Eqn == -fvc::grad(p1));

    fvOptions.correct(U1);
    K = 0.5*magSqr(U1);
}

G = mag(U1)*rhom;
