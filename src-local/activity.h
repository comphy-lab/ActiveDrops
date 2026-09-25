/**
# Emitted-product transport for active drops

The calling code defines `stracers` and supplies each tracer's dimensionless
diffusivity `D=1/Pe` and diffuse interfacial source coefficient
`A=AcNum/Pe`. The temporary VOF tracers transport complementary phase
contributions before reconstructing the concentration field. Diffusion is
weighted toward the selected phase through the tracer's `inverse` flag.

The scheme is based on Basilisk's `henry.h` pattern and the phase-change
transport of Farsoiya et al. (2021). Geometric PLIC fragments define the
discrete surface measure. Convergence of the coupled whole-field extension to
a sharp exterior-flux model is a separate verification problem.
*/

attribute {
  scalar phi1, phi2; // private
  double A;  // volumetric source coefficient AcNum/Pe
  double D; // dimensionless diffusivity 1/Pe in the diffusive phase
}

extern scalar * stracers;
scalar ActivityFlux[];

#include "diffusion.h"

/**
### activity_interface_area()

Returns the reconstructed interfacial area assigned to the current cell. PLIC
fragments are stored in mixed cells. If an interface lies exactly on a grid
face, no mixed cell exists; its face area is assigned once to the adjacent
pure exterior cell (`f=0`). In axisymmetry the returned area includes the
radial metric but omits the common factor $2\pi$.

This construction keeps the source out of pure drop cells and full solids.
It does not change the numerical extension used to store and transport the
concentration field.
*/
static inline double activity_interface_area (Point point, scalar phase)
{
  const double eps = 1e-6;
#if EMBED
  if (cs[] <= 0.)
    return 0.;
#endif

  if (phase[] > eps && phase[] < 1. - eps) {
    coord n = interface_normal(point, phase), p;
    double alpha = plane_alpha(phase[], n);
    double area = pow(Delta, dimension - 1)*
      plane_area_center(n, alpha, &p);
#if AXI
    area *= max(y + p.y*Delta, 0.);
#endif
    return area;
  }

  if (phase[] <= eps) {
    double area = 0.;
    foreach_dimension() {
      if (phase[-1] >= 1. - eps)
        area += fm.x[]*pow(Delta, dimension - 1);
      if (phase[1] >= 1. - eps)
        area += fm.x[1]*pow(Delta, dimension - 1);
    }
    return area;
  }

  return 0.;
}

/**
## Defaults

On trees we need to ensure conservation of the tracer when
refining/coarsening. */

event defaults (i = 0)
{
  for (scalar s in stracers) {
#if TREE
#if EMBED
s.refine = refine_embed_linear;
set_prolongation (s, refine_embed_linear);
#else
s.refine  = refine_bilinear;
#endif
s.restriction = restriction_volume_average;
s.gradient = p.gradient;
#endif // TREE
  }
}

/**
## Advection

To avoid numerical diffusion through the interface we use the [VOF
tracer transport scheme](/src/vof.h) for the temporary fields
$\phi_1$ and $\phi_2$, see section 3.2 of [Farsoiya et al.,
2021](#farsoiya2021). */

static scalar * phi_tracers = NULL;

event vof (i++)
{
  phi_tracers = f.tracers;
  for (scalar c in stracers) {
    scalar phi1 = new scalar, phi2 = new scalar;
    c.phi1 = phi1, c.phi2 = phi2;
    scalar_clone (phi1, c);
    scalar_clone (phi2, c);
    // Diffusion-side selection on c must not reverse the c*f VOF tracer.
    phi1.inverse = false;
    phi2.inverse = true;
    
    f.tracers = list_append (f.tracers, phi1);
    f.tracers = list_append (f.tracers, phi2);

    /**
    $\phi_1$ and $\phi_2$ are computed from $c$ as
    $$
    \phi_1 = c f
    $$
    $$
    \phi_2 = c (1-f)
    $$
    */
		  
    foreach() {
      double a = c[];
#if EMBED
      if (cs[] <= 0.)
        a = 0.;
#endif
      phi1[] = a*f[];
      phi2[] = a*(1. - f[]);
    }
  }
}

event tracer_diffusion (i++)
{
  free (f.tracers);
  f.tracers = phi_tracers;
  for (scalar c in stracers) {
    /**
    The advected concentration is computed from $\phi_1$ and $\phi_2$ as
    $$
    c = \phi_1 + \phi_2
    $$
    and these fields are then discarded. */
    
    scalar phi1 = c.phi1, phi2 = c.phi2;
    foreach() {
      c[] = phi1[] + phi2[];
#if EMBED
      if (cs[] <= 0.)
        c[] = 0.;
#endif
    }
    delete ({phi1, phi2});
    
    scalar volumic_metric[], dirichlet_source_term[];
    face vector diffusion_coef[];

    foreach() {
      double area = activity_interface_area(point, f);
      ActivityFlux[] = cm[] > 0. ?
        c.A*area/(cm[]*pow(Delta, dimension)) : 0.;
    }
  
    foreach() {
      volumic_metric[] = cm[];
      dirichlet_source_term[] = cm[]*ActivityFlux[];
    }
    foreach_face(){
      double ff = (f[] + f[-1])/2.;
      double wt = c.inverse ? (1.-ff) : ff;
      diffusion_coef.x[] = fm.x[]*c.D*wt;
    }
    diffusion (c, dt, D = diffusion_coef, r = dirichlet_source_term, theta = volumic_metric);
    }
}
