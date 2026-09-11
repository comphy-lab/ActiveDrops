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

For a circular planar drop with pure Marangoni forcing, choose the classical
mobility scales

$$
M_{2D}=\frac{R_0\gamma_C}{2(\mu_i+\mu_o)},\qquad
U_M=\frac{A_0M_{2D}}{D},\qquad C_A=\frac{A_0R_0}{D}.
$$

If $U_0=U_M$ and $C_*=C_A$, then $A_c=1$ and

$$
\Gamma_c=\frac{\gamma_C C_A}{\mu_oU_M}
=\frac{R_0\gamma_C}{\mu_oM_{2D}}=2(1+\lambda).
$$

The value four therefore follows from equal viscosities in this planar
normalization. It is independent of the spherical instability threshold
$Pe_M=4$. The positive slope follows from choosing an emitted product that
raises surface tension; mobility normalization determines its magnitude.
The associated unit outward-gradient condition is $-\partial_n c=1$,
so its dimensionless flux is $1/Pe$, not one.

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
Oh=\frac{\mu_o}{\sqrt{\rho_o\gamma_0R_0}}
=\sqrt{Ca/Re},\qquad Re=\frac{Ca}{Oh^2}.
$$

This identity follows directly from the reference scales. The momentum
coefficient is $Re$; a coefficient $4/Oh^2$ would not follow from this
definition of $Oh$ except at $Ca=4$. The planar mobility factor belongs in
the velocity and concentration normalization, not in an extra density factor.

With $\rho_r=1+(r_\rho-1)f$ and $\mu_r=1+(\lambda-1)f$, the momentum
equation used by the drivers is

$$
Re\,\rho_r(\partial_t\mathbf u+\mathbf u\cdot\nabla\mathbf u)
=-\nabla p+\nabla\cdot[\mu_r(\nabla\mathbf u+\nabla\mathbf u^T)]
+\mathbf f_\gamma,\qquad \nabla\cdot\mathbf u=0.
$$

Here $\mathbf f_\gamma$ is the interfacial surface-stress force. Both inertia
terms are retained. A comparison with an unsteady-Stokes approximation must
account for its omission of convective inertia, even at the same $Re$.

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
