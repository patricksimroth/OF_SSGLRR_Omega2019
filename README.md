# SSG/LRR-omega Reynolds-stress model for OpenFOAM

Implementation of the SSG/LRR-omega differential Reynolds-stress model with the 2019 length-scale
correction and simple diffusion (**SSGLRR-RSM-w2019-SD**) in OpenFOAM. The model solves six transport equations for the Reynolds stresses and one BSL omega equation.

The implementation follows the definition on the
[Turbulence Modeling Resource](https://tmbwg.github.io/turbmodels/rsm-ssglrr.html). The equations as implemented are collected in [equations.md](equations.md).

The model was developed with OpenFOAM v2512 (openfoam.com) and used with the compressible library on the 2D zero-pressure-gradient flat plate.

## Installation

To compile the model, execute:

```bash
./Allwmake
```

Once compiled, you should see `libSSGLRROmega2019CompressibleTurbulenceModel.so` in `$FOAM_USER_LIBBIN`.

To remove the compiled libraries and build files, execute:

```bash
./Allwclean
```

## Setting up a case

Add the library to `system/controlDict`:

```foam
libs    ("libSSGLRROmega2019CompressibleTurbulenceModel.so");
```

Select the model in `constant/turbulenceProperties`:

```foam
simulationType  RAS;

RAS
{
    RASModel        SSGLRROmega2019;
    turbulence      on;
    printCoeffs     on;
}
```

### Initial and boundary conditions

The model reads the fields `R` and `omega`. On a wall-resolved mesh, use:

```foam
// 0/R
wall
{
    type            fixedValue;
    value           uniform (0 0 0 0 0 0);
}

// 0/omega
wall
{
    type            omegaWallFunction;
    blending        stepwise;
    beta1           0.0075;
    value           $internalField;
}
```

The Turbulence Modeling Resource prescribes `omega_wall = 10*6*nu/(beta_1*d^2)` with `beta_1 = 0.075`.
`omegaWallFunction` evaluates `6*nu/(beta1*d^2)` below the laminar sublayer limit, so `beta1 0.0075`
reproduces that value.

### Discretisation schemes

The model needs the following entries in `system/fvSchemes`, in addition to the schemes of the
flow solver:

```foam
divSchemes
{
    div(phi,R)                                      bounded Gauss limitedLinear 1;
    div(phi,omega)                                  bounded Gauss limitedLinear 1;
    div((rho*R))                                    Gauss linear;
    div(((rho*(thermo:mu|rho))*dev2(T(grad(U)))))   Gauss linear;
    div(devRhoReff)                                 Gauss linear;
}
```

`bounded` is only needed for RANS. For URANS, it can be left out.

## Model coefficients

All coefficients are optional. The defaults below are the published values, and they can be changed in the
`SSGLRROmega2019Coeffs` dictionary:

```foam
RAS
{
    RASModel        SSGLRROmega2019;

    SSGLRROmega2019Coeffs
    {
        Cmu             0.09;
        kappa           0.41;

        // Length-scale correction
        lscA            31;
        lscChiT         1;

        // Omega equation
        alphaOmegaLRR   0.5556;
        alphaOmegaSSG   0.44;
        betaOmegaLRR    0.075;
        betaOmegaSSG    0.0828;
        sigmaOmegaLRR   0.5;
        sigmaOmegaSSG   0.856;
        sigmaDSSG       1.712;

        // Pressure-strain correlation
        C1LRR           1.8;
        C1SSG           1.7;
        C1sSSG          0.9;
        C2LRR           0.52;
        C2SSG           1.05;
        C3              0.8;
        C3sSSG          0.65;
        C4SSG           0.625;
        C5SSG           0.2;

        // Simple diffusion
        DLRR            0.5;
        DSSG            0.22;
    }
}
```

## Differences from the published model

- The production in the omega equation uses `G = 0.5*mag(tr(P))` instead of `0.5*tr(P)`, which keeps it
  non-negative.
- After the normal stresses have been bounded, the shear stresses are clipped to the realizable range
  `|Rij| <= sqrt(Rii*Rjj)`. Without this limit, cells with clipped normal stresses can diverge on
  unstructured meshes.

## References

- Turbulence Modeling Resource, *SSG/LRR-omega full Reynolds stress model*,
  [https://tmbwg.github.io/turbmodels/rsm-ssglrr.html](https://tmbwg.github.io/turbmodels/rsm-ssglrr.html)
- B. Eisfeld, C. Rumsey, and V. Togiti, *Verification and Validation of a Second-Moment-Closure Model*,
  AIAA Journal, Vol. 54, no. 5, pp. 1524-1541, 2016, [https://doi.org/10.2514/1.J054718](https://doi.org/10.2514/1.J054718)
- R.-D. Cécora, R. Radespiel, B. Eisfeld, and A. Probst, *Differential Reynolds-Stress Modeling for
  Aeronautics*, AIAA Journal, Vol. 53, no. 3, pp. 739-755, 2015,
  [https://doi.org/10.2514/1.J053250](https://doi.org/10.2514/1.J053250)
- B. Eisfeld and C. L. Rumsey, *Length-Scale Correction for Reynolds-Stress Modeling*, AIAA Journal,
  Vol. 58, no. 4, pp. 1518-1528, 2020, [https://doi.org/10.2514/1.J058858](https://doi.org/10.2514/1.J058858)
- The repository structure follows [OF_SA_SALSA](https://github.com/JanisGeise/OF_SA_SALSA).

## License

This project is licensed under the GNU General Public License v3.0, see [LICENSE](LICENSE).
