/**
# Dimensionless model for solubilizing active drops

The model uses the initial drop radius $R_0$, a chosen velocity $U_0$ and a
chosen emitted-species concentration $C_*$ as independent scales. Viscosity
is scaled by $\mu_o$, pressure and stress by $\mu_oU_0/R_0$, and surface
tension by $\mu_oU_0$. The input groups are

$$
Re = \frac{\rho_o U_0R_0}{\mu_o},\qquad
Ca = \frac{\mu_oU_0}{\gamma_0},\qquad
Pe = \frac{U_0R_0}{D},
$$

$$
\Gamma_c = \frac{\gamma_C C_*}{\mu_oU_0},\qquad
A_c = \frac{A_0R_0}{DC_*},\qquad
\lambda = \frac{\mu_i}{\mu_o},\qquad
r_\rho = \frac{\rho_i}{\rho_o}.
$$

Here $A_0$ is the dimensional outward flux of a coarse-grained product that
increases surface tension, and $\gamma_C = \partial\gamma/\partial C > 0$.
The implemented dimensionless equations therefore use

$$
D_c = Pe^{-1},\qquad
S_c = \frac{A_c}{Pe}\,\delta_{\Gamma,h},\qquad
\gamma = Ca^{-1} + \Gamma_c c.
$$

This is a prescribed-flux solubilization model.  It conserves the drop phase
approximately through CLSVOF and does not resolve micelle formation, droplet
mass loss or a finite fuel inventory.

The discrete surface measure $\delta_{\Gamma,h}$ is the geometric PLIC
fragment area divided by the metric cell volume.  Mixed cells carry their
reconstructed fragment; an exactly face-aligned interface is assigned to its
pure exterior neighbour.  Pure drop cells and full solids receive no source.
The concentration remains a whole-field numerical extension, with diffusion
weighted toward the outer phase.

## Mobility diagnostics

The material coupling $\Gamma_c$ is independent of geometry.  Geometry enters
only when comparing the chosen velocity scale with the classical mobility
scale.  For a circular planar drop and a spherical drop, respectively,

$$
G_{2D}=2(1+\lambda),\qquad G_{3D}=2+3\lambda,
$$

$$
\chi = \frac{U_M}{U_0} = \frac{A_c\Gamma_c}{G},\qquad
Pe_M = \chi Pe.
$$

Thus $\Gamma_c=4$ at $A_c=\lambda=1$ normalizes the planar mobility, while a
spherical mobility normalization would use $5$.  The simulations retain the
same default material value $\Gamma_c=4$ in all geometries and report $\chi$
and $Pe_M$ rather than silently changing the material in the pipe.

The reference-tension Ohnesorge number is derived, never independently specified:

$$
Oh = \sqrt{Ca/Re}.
$$

## Evidence boundary

The scales follow Michelin (2023), equations (2)--(5), and the planar mobility
used by Li & Koch (2022), DOI 10.1017/jfm.2022.891.  The spherical and planar
mobility formulae are comparison diagnostics.  They do not verify the diffuse
source, predict a confined onset or turn a finite-time `MOVED` classification
into a critical Péclet number.
*/

#ifndef ACTIVE_DROP_MODEL_H
#define ACTIVE_DROP_MODEL_H

#include <math.h>

static inline double active_drop_ohnesorge (double reynolds,
                                             double capillary)
{
  return sqrt(capillary/reynolds);
}

static inline double active_drop_mobility_denominator
  (double viscosity_ratio, int spherical)
{
  return spherical ? 2. + 3.*viscosity_ratio :
    2.*(1. + viscosity_ratio);
}

static inline double active_drop_velocity_scale_ratio
  (double activity, double gamma_slope, double viscosity_ratio,
   int spherical)
{
  return activity*gamma_slope/
    active_drop_mobility_denominator(viscosity_ratio, spherical);
}

#endif
