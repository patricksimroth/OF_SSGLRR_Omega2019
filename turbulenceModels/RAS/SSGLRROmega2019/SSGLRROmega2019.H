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

Class
    Foam::RASModels::SSGLRROmega2019

Group
    grpRASTurbulence

Description
    SSG/LRR-omega differential Reynolds-stress turbulence model with the 2019
    length-scale correction and simple diffusion (SSGLRR-RSM-w2019-SD) for
    compressible flows.

    The coefficients blend from the LRR set near the wall to the SSG set in
    the outer flow with Menter's BSL function F1, phi = F1*phiLRR
    + (1 - F1)*phiSSG.  The length scale is provided by the BSL omega equation.

    Required fields
    \verbatim
        R         | Reynolds stress tensor [m2/s2]
        omega     | Specific dissipation rate [1/s]
    \endverbatim

    References:
    \verbatim
        Model definition:
            Turbulence Modeling Resource (2026).
            SSG/LRR-omega full Reynolds stress model.
            https://tmbwg.github.io/turbmodels/rsm-ssglrr.html
            (Retrieved:18-11-2026).

        SSG/LRR-omega model (w2012):
            Eisfeld, B., Rumsey, C., & Togiti, V. (2016).
            Verification and validation of a second-moment-closure model.
            AIAA Journal, 54(5), 1524-1541.
            DOI:10.2514/1.J054718

            Cecora, R.-D., Radespiel, R., Eisfeld, B., & Probst, A. (2015).
            Differential Reynolds-stress modeling for aeronautics.
            AIAA Journal, 53(3), 739-755.
            DOI:10.2514/1.J053250

        Length-scale correction (w2019):
            Eisfeld, B., & Rumsey, C. L. (2020).
            Length-scale correction for Reynolds-stress modeling.
            AIAA Journal, 58(4), 1518-1528.
            DOI:10.2514/1.J058858
    \endverbatim

Usage
    Example by using \c constant/turbulenceProperties:
    \verbatim
    RAS
    {
        // Mandatory entries
        RASModel        SSGLRROmega2019;

        // Optional entries
        turbulence      on;
        printCoeffs     on;

        // Default coefficients
        SSGLRROmega2019Coeffs
        {
            Cmu             0.09;
            kappa           0.41;

            lscA            31;
            lscChiT         1;

            alphaOmegaLRR   0.5556;
            alphaOmegaSSG   0.44;
            betaOmegaLRR    0.075;
            betaOmegaSSG    0.0828;
            sigmaOmegaLRR   0.5;
            sigmaOmegaSSG   0.856;
            sigmaDSSG       1.712;

            C1LRR           1.8;
            C1SSG           1.7;
            C1sSSG          0.9;
            C2LRR           0.52;
            C2SSG           1.05;
            C3              0.8;
            C3sSSG          0.65;
            C4SSG           0.625;
            C5SSG           0.2;

            DLRR            0.5;
            DSSG            0.22;
        }
    }
    \endverbatim

Note
  - C1s, C2, C3s and sigmaD are zero in the LRR set.  C4 and C5 of the LRR set
    are derived from \c C2LRR.  The simple-diffusion coefficient is
    D = F1*DLRR*Cmu + (1 - F1)*(2/3)*DSSG.
  - The production in the omega equation is limited to G = 0.5*mag(tr(P)),
    which keeps it non-negative.
  - The shear stresses are clipped to the realizable range
    |Rij| <= sqrt(Rii*Rjj) after the normal stresses have been bounded.

SourceFiles
    SSGLRROmega2019.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_SSGLRROmega2019_H
#define Foam_SSGLRROmega2019_H

#include "RASModel.H"
#include "ReynoldsStress.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                       Class SSGLRROmega2019 Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class SSGLRROmega2019
:
    public ReynoldsStress<RASModel<BasicTurbulenceModel>>
{
    // Private Member Functions

        //- No copy construct
        SSGLRROmega2019(const SSGLRROmega2019&) = delete;

        //- No copy assignment
        void operator=(const SSGLRROmega2019&) = delete;


protected:

    // Protected Data

        // Model constants

            //- Eddy-viscosity coefficient
            dimensionedScalar Cmu_;

            //- von Karman constant
            dimensionedScalar kappa_;


            // Length-scale correction

                //- Slope of the correction function
                dimensionedScalar lscA_;

                //- Threshold of the correction function
                dimensionedScalar lscChiT_;


            // Omega equation

                //- Production coefficient of the LRR set
                dimensionedScalar alphaOmegaLRR_;

                //- Production coefficient of the SSG set
                dimensionedScalar alphaOmegaSSG_;

                //- Destruction coefficient of the LRR set
                dimensionedScalar betaOmegaLRR_;

                //- Destruction coefficient of the SSG set
                dimensionedScalar betaOmegaSSG_;

                //- Diffusion coefficient of the LRR set
                dimensionedScalar sigmaOmegaLRR_;

                //- Diffusion coefficient of the SSG set
                dimensionedScalar sigmaOmegaSSG_;

                //- Cross-diffusion coefficient of the SSG set
                dimensionedScalar sigmaDSSG_;


            // Pressure-strain correlation

                //- Slow-term coefficient C1 of the LRR set
                dimensionedScalar C1LRR_;

                //- Slow-term coefficient C1 of the SSG set
                dimensionedScalar C1SSG_;

                //- Production-dependent coefficient C1* of the SSG set
                dimensionedScalar C1sSSG_;

                //- LRR coefficient from which C4 and C5 of the LRR set follow
                dimensionedScalar C2LRR_;

                //- Quadratic-term coefficient C2 of the SSG set
                dimensionedScalar C2SSG_;

                //- Rapid-term coefficient C3 of both sets
                dimensionedScalar C3_;

                //- Anisotropy-dependent coefficient C3* of the SSG set
                dimensionedScalar C3sSSG_;

                //- Strain-term coefficient C4 of the SSG set
                dimensionedScalar C4SSG_;

                //- Rotation-term coefficient C5 of the SSG set
                dimensionedScalar C5SSG_;


            // Simple diffusion

                //- Diffusion coefficient of the LRR set, multiplies Cmu
                dimensionedScalar DLRR_;

                //- Diffusion coefficient of the SSG set, multiplies 2/3
                dimensionedScalar DSSG_;


        // Fields

            //- Turbulent kinetic energy [m2/s2]
            volScalarField k_;

            //- Specific dissipation rate [1/s]
            volScalarField omega_;

            //- Wall distance [m]
            //  Note: different to wall distance in parent RASModel
            //  which is for near-wall cells only
            const volScalarField& y_;


    // Protected Member Functions

        //- Return the limited product grad(k) & grad(omega)
        tmp<volScalarField> gradKGradOmega() const;

        //- Return the blending function between the LRR and the SSG set
        tmp<volScalarField> F1() const;

        //- Return the length-scale correction function
        tmp<volScalarField> FLSC() const;

        //- Return the coefficient blended between the LRR and the SSG set
        tmp<volScalarField> blend
        (
            const volScalarField& F1,
            const dimensionedScalar& psiLRR,
            const dimensionedScalar& psiSSG
        ) const;

        //- Clip the shear stresses to |Rij| <= sqrt(Rii*Rjj)
        void boundShearStress(volSymmTensorField& R) const;

        //- Update the eddy viscosity
        virtual void correctNut();


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("SSGLRROmega2019");


    // Constructors

        //- Construct from components
        SSGLRROmega2019
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~SSGLRROmega2019() = default;


    // Member Functions

        //- Read model coefficients if they have changed
        virtual bool read();

        //- Return the turbulent kinetic energy
        virtual tmp<volScalarField> k() const;

        //- Return the turbulent kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const;

        //- Solve the omega and the Reynolds-stress equations
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "SSGLRROmega2019.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
