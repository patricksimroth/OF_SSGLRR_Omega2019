# Source-Code Consistent Equations

This document summarizes the equations implemented by `SSGLRROmega2019`. They are written in the form used
in the code, with the kinematic Reynolds stress $R_{ij}$ and the turbulent kinetic energy $k = R_{ii}/2$.
OpenFOAM's phase-fraction factor and the source terms from `fvOptions` are omitted. For incompressible usage,
the density $\rho$ is effectively one.

## Blending

Coefficients blend between the LRR set near the wall and the SSG set in the outer flow:

$$
\phi = F_1 \phi^{LRR} + (1 - F_1) \phi^{SSG},
$$

with

$$
F_1 = \tanh(\zeta^4),
$$

$$
\zeta = \min\left[\max\left(\frac{\sqrt{k}}{C_\mu \omega d}, \frac{500 \nu}{\omega d^2}\right),
\frac{4 \sigma_\omega^{SSG} k}{CD_{k\omega} d^2}\right],
$$

$$
CD_{k\omega} = \max\left[\frac{\sigma_d^{SSG}}{\omega} \max\left(\frac{\partial k}{\partial x_j}
\frac{\partial \omega}{\partial x_j}, 0\right), 10^{-15}\right],
$$

where $d$ is the distance to the nearest wall.

## Reynolds-stress equation

$$
\frac{\partial (\rho R_{ij})}{\partial t} + \frac{\partial (\rho U_k R_{ij})}{\partial x_k} =
\rho P_{ij} + \rho \Pi_{ij} - \frac{2}{3} \rho \varepsilon \delta_{ij} +
\frac{\partial}{\partial x_k}\left[\rho \left(\nu + \frac{D}{C_\mu} \frac{k}{\omega}\right)
\frac{\partial R_{ij}}{\partial x_k}\right],
$$

with the production and the dissipation rate

$$
P_{ij} = -R_{ik} \frac{\partial U_j}{\partial x_k} - R_{jk} \frac{\partial U_i}{\partial x_k},
\qquad
\varepsilon = C_\mu k \omega.
$$

The pressure-strain correlation is

$$
\begin{aligned}
\Pi_{ij} = {} & -\left(C_1 \varepsilon + C_1^* G\right) a_{ij}
+ C_2 \varepsilon \left(a_{ik} a_{kj} - \frac{1}{3} a_{kl} a_{kl} \delta_{ij}\right) \\
& + \left(C_3 - C_3^* \sqrt{a_{kl} a_{kl}}\right) k S_{ij}^* \\
& + C_4 k \left(a_{ik} S_{jk} + a_{jk} S_{ik} - \frac{2}{3} a_{kl} S_{kl} \delta_{ij}\right)
+ C_5 k \left(a_{ik} W_{jk} + a_{jk} W_{ik}\right),
\end{aligned}
$$

with

$$
a_{ij} = \frac{R_{ij}}{k} - \frac{2}{3} \delta_{ij},
\qquad
G = \frac{1}{2} \left|P_{kk}\right|,
$$

$$
S_{ij} = \frac{1}{2}\left(\frac{\partial U_i}{\partial x_j} + \frac{\partial U_j}{\partial x_i}\right),
\qquad
S_{ij}^* = S_{ij} - \frac{1}{3} S_{kk} \delta_{ij},
\qquad
W_{ij} = \frac{1}{2}\left(\frac{\partial U_i}{\partial x_j} - \frac{\partial U_j}{\partial x_i}\right).
$$

The simple-diffusion coefficient is

$$
D = F_1 D^{LRR} C_\mu + (1 - F_1) \frac{2}{3} D^{SSG}.
$$

## Omega equation

$$
\frac{\partial (\rho \omega)}{\partial t} + \frac{\partial (\rho U_k \omega)}{\partial x_k} =
\rho \alpha_\omega \frac{\omega}{k} G - \left(1 - F^{LSC}\right) \rho \beta_\omega \omega^2 +
\frac{\partial}{\partial x_k}\left[\rho \left(\nu + \sigma_\omega \frac{k}{\omega}\right)
\frac{\partial \omega}{\partial x_k}\right] +
\rho \frac{\sigma_d}{\omega} \max\left(\frac{\partial k}{\partial x_j}
\frac{\partial \omega}{\partial x_j}, 0\right).
$$

The length-scale correction is

$$
F^{LSC} = \frac{1}{2}\left\{1 + \tanh\left[A \left(\chi - \chi_T\right)\right]\right\},
$$

$$
\chi = \max\left[\left(\frac{L}{\ell} - 1\right) \left(\frac{L}{\ell}\right)^2, 0\right],
\qquad
\frac{L}{\ell} = \frac{\sqrt{k}}{C_\mu^{1/4} \kappa \omega d}.
$$

## Eddy viscosity

The eddy viscosity is

$$
\nu_t = \frac{k}{\omega}.
$$

## Limiters

After each solution of the Reynolds-stress equation, the normal stresses are bounded with OpenFOAM's
`boundNormalStress`, and the shear stresses are clipped to

$$
\left|R_{ij}\right| \le \sqrt{R_{ii} R_{jj}} \qquad (i \ne j, \text{ no summation}).
$$

## Coefficients

| Coefficient | Dictionary entries | LRR set | SSG set |
|---|---|---|---|
| $\alpha_\omega$ | `alphaOmegaLRR`, `alphaOmegaSSG` | 0.5556 | 0.44 |
| $\beta_\omega$ | `betaOmegaLRR`, `betaOmegaSSG` | 0.075 | 0.0828 |
| $\sigma_\omega$ | `sigmaOmegaLRR`, `sigmaOmegaSSG` | 0.5 | 0.856 |
| $\sigma_d$ | `sigmaDSSG` | 0 | 1.712 |
| $C_1$ | `C1LRR`, `C1SSG` | 1.8 | 1.7 |
| $C_1^*$ | `C1sSSG` | 0 | 0.9 |
| $C_2$ | `C2SSG` | 0 | 1.05 |
| $C_3$ | `C3` | 0.8 | 0.8 |
| $C_3^*$ | `C3sSSG` | 0 | 0.65 |
| $C_4$ | `C2LRR`, `C4SSG` | $(9 C_2^{LRR} + 6)/11$ | 0.625 |
| $C_5$ | `C2LRR`, `C5SSG` | $(-7 C_2^{LRR} + 10)/11$ | 0.2 |
| $D$ factors | `DLRR`, `DSSG` | 0.5 | 0.22 |

The remaining constants are $C_\mu = 0.09$ (`Cmu`), $\kappa = 0.41$ (`kappa`), $C_2^{LRR} = 0.52$ (`C2LRR`),
$A = 31$ (`lscA`), and $\chi_T = 1$ (`lscChiT`).
