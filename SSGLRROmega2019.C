/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
    Copyright (C) 2026 Patrick Simroth
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

#include "SSGLRROmega2019.H"
#include "wallDist.H"
#include "fvOptions.H"
#include "bound.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField>
SSGLRROmega2019<BasicTurbulenceModel>::gradKGradOmega() const
{
    return max
    (
        fvc::grad(k_) & fvc::grad(omega_),
        dimensionedScalar(dimless/pow3(dimTime), 0)
    );
}


template<class BasicTurbulenceModel>
tmp<volScalarField> SSGLRROmega2019<BasicTurbulenceModel>::F1() const
{
    const volScalarField CDkOmega
    (
        max
        (
            sigmaDSSG_*gradKGradOmega()/omega_,
            dimensionedScalar(sqr(dimless/dimTime), SMALL)
        )
    );

    const volScalarField zeta
    (
        min
        (
            max
            (
                sqrt(k_)/(Cmu_*omega_*y_),
                500*this->nu()/(omega_*sqr(y_))
            ),
            4*sigmaOmegaSSG_*k_/(CDkOmega*sqr(y_))
        )
    );

    return tanh(pow4(zeta));
}


template<class BasicTurbulenceModel>
tmp<volScalarField> SSGLRROmega2019<BasicTurbulenceModel>::FLSC() const
{
    // Ratio of the turbulent length scale to the log-layer length scale
    const volScalarField ratio
    (
        sqrt(k_)/(pow(Cmu_, 0.25)*kappa_*omega_*y_)
    );

    const volScalarField chi(max((ratio - scalar(1))*sqr(ratio), scalar(0)));

    return 0.5*(scalar(1) + tanh(lscA_*(chi - lscChiT_)));
}


template<class BasicTurbulenceModel>
tmp<volScalarField> SSGLRROmega2019<BasicTurbulenceModel>::blend
(
    const volScalarField& F1,
    const dimensionedScalar& psiLRR,
    const dimensionedScalar& psiSSG
) const
{
    return psiSSG + F1*(psiLRR - psiSSG);
}


template<class BasicTurbulenceModel>
void SSGLRROmega2019<BasicTurbulenceModel>::boundShearStress
(
    volSymmTensorField& R
) const
{
    symmTensorField& Rcells = R.primitiveFieldRef();

    forAll(Rcells, celli)
    {
        symmTensor& Rc = Rcells[celli];

        const scalar maxXY = sqrt(Rc.xx()*Rc.yy());
        const scalar maxXZ = sqrt(Rc.xx()*Rc.zz());
        const scalar maxYZ = sqrt(Rc.yy()*Rc.zz());

        Rc.xy() = max(min(Rc.xy(), maxXY), -maxXY);
        Rc.xz() = max(min(Rc.xz(), maxXZ), -maxXZ);
        Rc.yz() = max(min(Rc.yz(), maxYZ), -maxYZ);
    }

    R.correctBoundaryConditions();
}


template<class BasicTurbulenceModel>
void SSGLRROmega2019<BasicTurbulenceModel>::correctNut()
{
    this->nut_ = k_/omega_;
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
SSGLRROmega2019<BasicTurbulenceModel>::SSGLRROmega2019
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    ReynoldsStress<RASModel<BasicTurbulenceModel>>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    Cmu_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cmu",
            this->coeffDict_,
            0.09
        )
    ),
    kappa_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "kappa",
            this->coeffDict_,
            0.41
        )
    ),

    lscA_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "lscA",
            this->coeffDict_,
            31
        )
    ),
    lscChiT_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "lscChiT",
            this->coeffDict_,
            1
        )
    ),

    alphaOmegaLRR_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "alphaOmegaLRR",
            this->coeffDict_,
            0.5556
        )
    ),
    alphaOmegaSSG_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "alphaOmegaSSG",
            this->coeffDict_,
            0.44
        )
    ),
    betaOmegaLRR_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "betaOmegaLRR",
            this->coeffDict_,
            0.075
        )
    ),
    betaOmegaSSG_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "betaOmegaSSG",
            this->coeffDict_,
            0.0828
        )
    ),
    sigmaOmegaLRR_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "sigmaOmegaLRR",
            this->coeffDict_,
            0.5
        )
    ),
    sigmaOmegaSSG_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "sigmaOmegaSSG",
            this->coeffDict_,
            0.856
        )
    ),
    sigmaDSSG_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "sigmaDSSG",
            this->coeffDict_,
            1.712
        )
    ),

    C1LRR_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "C1LRR",
            this->coeffDict_,
            1.8
        )
    ),
    C1SSG_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "C1SSG",
            this->coeffDict_,
            1.7
        )
    ),
    C1sSSG_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "C1sSSG",
            this->coeffDict_,
            0.9
        )
    ),
    C2LRR_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "C2LRR",
            this->coeffDict_,
            0.52
        )
    ),
    C2SSG_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "C2SSG",
            this->coeffDict_,
            1.05
        )
    ),
    C3_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "C3",
            this->coeffDict_,
            0.8
        )
    ),
    C3sSSG_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "C3sSSG",
            this->coeffDict_,
            0.65
        )
    ),
    C4SSG_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "C4SSG",
            this->coeffDict_,
            0.625
        )
    ),
    C5SSG_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "C5SSG",
            this->coeffDict_,
            0.2
        )
    ),

    DLRR_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "DLRR",
            this->coeffDict_,
            0.5
        )
    ),
    DSSG_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "DSSG",
            this->coeffDict_,
            0.22
        )
    ),

    k_
    (
        IOobject
        (
            "k",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        0.5*tr(this->R_)
    ),
    omega_
    (
        IOobject
        (
            "omega",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    y_(wallDist::New(this->mesh_).y())
{
    if (type == typeName)
    {
        this->printCoeffs(type);

        this->boundNormalStress(this->R_);
        bound(omega_, this->omegaMin_);
        k_ = 0.5*tr(this->R_);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool SSGLRROmega2019<BasicTurbulenceModel>::read()
{
    if (ReynoldsStress<RASModel<BasicTurbulenceModel>>::read())
    {
        Cmu_.readIfPresent(this->coeffDict());
        kappa_.readIfPresent(this->coeffDict());

        lscA_.readIfPresent(this->coeffDict());
        lscChiT_.readIfPresent(this->coeffDict());

        alphaOmegaLRR_.readIfPresent(this->coeffDict());
        alphaOmegaSSG_.readIfPresent(this->coeffDict());
        betaOmegaLRR_.readIfPresent(this->coeffDict());
        betaOmegaSSG_.readIfPresent(this->coeffDict());
        sigmaOmegaLRR_.readIfPresent(this->coeffDict());
        sigmaOmegaSSG_.readIfPresent(this->coeffDict());
        sigmaDSSG_.readIfPresent(this->coeffDict());

        C1LRR_.readIfPresent(this->coeffDict());
        C1SSG_.readIfPresent(this->coeffDict());
        C1sSSG_.readIfPresent(this->coeffDict());
        C2LRR_.readIfPresent(this->coeffDict());
        C2SSG_.readIfPresent(this->coeffDict());
        C3_.readIfPresent(this->coeffDict());
        C3sSSG_.readIfPresent(this->coeffDict());
        C4SSG_.readIfPresent(this->coeffDict());
        C5SSG_.readIfPresent(this->coeffDict());

        DLRR_.readIfPresent(this->coeffDict());
        DSSG_.readIfPresent(this->coeffDict());

        return true;
    }

    return false;
}


template<class BasicTurbulenceModel>
tmp<volScalarField> SSGLRROmega2019<BasicTurbulenceModel>::k() const
{
    return k_;
}


template<class BasicTurbulenceModel>
tmp<volScalarField> SSGLRROmega2019<BasicTurbulenceModel>::epsilon() const
{
    return Cmu_*k_*omega_;
}


template<class BasicTurbulenceModel>
void SSGLRROmega2019<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    volSymmTensorField& R = this->R_;
    fv::options& fvOptions = fv::options::New(this->mesh_);

    ReynoldsStress<RASModel<BasicTurbulenceModel>>::correct();

    const volTensorField gradU(fvc::grad(U));

    // Reynolds-stress production
    const volSymmTensorField P(-twoSymm(R & gradU));
    const volScalarField G(this->GName(), 0.5*mag(tr(P)));

    const volScalarField F1(this->F1());

    // Update omega and G at the wall
    omega_.boundaryFieldRef().updateCoeffs();
    omega_.boundaryFieldRef().template evaluateCoupled<coupledFvPatch>();

    // Blended coefficients of the omega equation
    const volScalarField alphaOmega(blend(F1, alphaOmegaLRR_, alphaOmegaSSG_));
    const volScalarField betaOmega(blend(F1, betaOmegaLRR_, betaOmegaSSG_));
    const volScalarField sigmaOmega(blend(F1, sigmaOmegaLRR_, sigmaOmegaSSG_));
    const volScalarField sigmaD(sigmaDSSG_*(scalar(1) - F1));

    const volScalarField crossDiffusion(sigmaD/omega_*gradKGradOmega());
    const volScalarField FLSC(this->FLSC());

    // Specific dissipation rate equation
    tmp<fvScalarMatrix> omegaEqn
    (
        fvm::ddt(alpha, rho, omega_)
      + fvm::div(alphaRhoPhi, omega_)
      - fvm::laplacian(alpha*rho*(this->nu() + sigmaOmega*k_/omega_), omega_)
     ==
        alpha*rho*alphaOmega*omega_/k_*G
      - fvm::Sp(alpha*rho*(scalar(1) - FLSC)*betaOmega*omega_, omega_)
      + alpha*rho*crossDiffusion
      + fvOptions(alpha, rho, omega_)
    );

    omegaEqn.ref().relax();
    fvOptions.constrain(omegaEqn.ref());
    omegaEqn.ref().boundaryManipulate(omega_.boundaryFieldRef());
    solve(omegaEqn);
    fvOptions.correct(omega_);
    bound(omega_, this->omegaMin_);

    const volScalarField epsilon(Cmu_*k_*omega_);

    // Anisotropy, strain-rate and rotation-rate tensors
    const volSymmTensorField a(dev(R)/k_);
    const volSymmTensorField S(symm(gradU));
    const volTensorField W(skew(gradU));

    // C4 and C5 of the LRR set
    const dimensionedScalar C4LRR = (9*C2LRR_ + 6.0)/11.0;
    const dimensionedScalar C5LRR = (-7*C2LRR_ + 10.0)/11.0;

    // Blended coefficients of the Reynolds-stress equation
    const volScalarField C1(blend(F1, C1LRR_, C1SSG_));
    const volScalarField C1s(C1sSSG_*(scalar(1) - F1));
    const volScalarField C2(C2SSG_*(scalar(1) - F1));
    const volScalarField C3s(C3sSSG_*(scalar(1) - F1));
    const volScalarField C4(blend(F1, C4LRR, C4SSG_));
    const volScalarField C5(blend(F1, C5LRR, C5SSG_));
    const volScalarField D
    (
        DLRR_*Cmu_*F1 + (2.0/3.0)*DSSG_*(scalar(1) - F1)
    );

    // Slow part of the pressure-strain correlation, -A*a
    const volScalarField A(C1*epsilon + C1s*G);

    // Reynolds-stress equation
    tmp<fvSymmTensorMatrix> REqn
    (
        fvm::ddt(alpha, rho, R)
      + fvm::div(alphaRhoPhi, R)
      - fvm::laplacian(alpha*rho*(this->nu() + D/Cmu_*k_/omega_), R)
      + fvm::Sp(alpha*rho*A/k_, R)
     ==
        alpha*rho
       *(
            P
          + (2.0/3.0)*A*I
          - (2.0/3.0)*epsilon*I
          + C2*epsilon*dev(innerSqr(a))
          + k_
           *(
                (C3_ - C3s*mag(a))*dev(S)
              + C4*devTwoSymm(a & S)
              + C5*twoSymm(a & W)
            )
        )
      + fvOptions(alpha, rho, R)
    );

    REqn.ref().relax();
    fvOptions.constrain(REqn.ref());
    solve(REqn);
    fvOptions.correct(R);

    this->boundNormalStress(R);
    boundShearStress(R);

    k_ = 0.5*tr(R);
    bound(k_, this->kMin_);

    correctNut();

    this->correctWallShearStress(R);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
