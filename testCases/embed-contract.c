/**
# Embedded geometry and impermeable scalar contracts

Small software checks for straight planar and axisymmetric walls after
repeated tree refinement and coarsening. Geometric integrals are evaluated
independently of the construction routine. A constant scalar must remain
constant under diffusion with homogeneous Neumann conditions. These checks
do not establish active-drop dynamics or a contact-line model.
*/

#include "grid/quadtree.h"
#include "embed.h"
#if TEST_PIPE
#include "axi.h"
#endif
#include "run.h"
#include "diffusion.h"
#include "embed-channel-geometry.h"

scalar concentration[];
scalar absorbing[];
concentration[embed] = neumann(0.);
concentration[left] = neumann(0.);
concentration[right] = neumann(0.);
concentration[top] = neumann(0.);
concentration[bottom] = neumann(0.);
absorbing[embed] = dirichlet(0.);
absorbing[left] = neumann(0.);
absorbing[right] = neumann(0.);
absorbing[top] = neumann(0.);
absorbing[bottom] = neumann(0.);

static const double wall = 0.731;
static int failed = 0;

/**
The antiderivative integrates a radial power from the lower fluid boundary
to a clamped coordinate. This also evaluates exact face integrals.
*/
static double primitive (double position, int radial)
{
  double lower = TEST_PIPE ? 0. : -wall;
  double z = fmin(fmax(position, lower), wall);
  return radial ? (z*z - lower*lower)/2. : z - lower;
}

static int check_geometry (int cycle)
{
  double volume = 0., cell_error = 0., face_error = 0.;
  int cut_cells = 0;
  foreach (reduction(+:volume) reduction(max:cell_error)
           reduction(+:cut_cells)) {
    double lo = y - Delta/2., hi = y + Delta/2.;
    double exact_cs = (primitive(hi, 0) - primitive(lo, 0))/Delta;
    double exact_cm = (primitive(hi, TEST_PIPE) -
                       primitive(lo, TEST_PIPE))/Delta;
    cell_error = max(cell_error, fabs(cs[] - exact_cs));
    cell_error = max(cell_error, fabs(cm[] - exact_cm));
    volume += cm[]*sq(Delta);
    if (cs[] > 0. && cs[] < 1.) cut_cells++;
  }
  foreach_face(x, reduction(max:face_error)) {
    double lo = y - Delta/2., hi = y + Delta/2.;
    double fraction = (primitive(hi, 0) - primitive(lo, 0))/Delta;
    double metric = (primitive(hi, TEST_PIPE) -
                     primitive(lo, TEST_PIPE))/Delta;
    face_error = max(face_error, fabs(fs.x[] - fraction));
    face_error = max(face_error, fabs(fm.x[] - metric));
  }
  foreach_face(y, reduction(max:face_error)) {
    double fraction = y < wall && (TEST_PIPE || y > -wall) ? 1. : 0.;
    double metric = TEST_PIPE ? fraction*max(y, 1e-20) : fraction;
    face_error = max(face_error, fabs(fs.y[] - fraction));
    face_error = max(face_error, fabs(fm.y[] - metric));
  }
  double exact_volume = L0*(TEST_PIPE ? wall*wall/2. : 2.*wall);
  int bad = !isfinite(volume) || fabs(volume - exact_volume) > 1e-11 ||
    cell_error > 1e-11 || face_error > 1e-11 || cut_cells == 0;
  fprintf(stderr, "geometry pipe=%d cycle=%d volume_error=%g cell_error=%g "
          "face_error=%g cuts=%d\n", TEST_PIPE, cycle,
          fabs(volume - exact_volume), cell_error, face_error, cut_cells);
  return bad;
}

static int check_diffusion (int cycle)
{
  scalar capacity[];
  face vector diffusivity[];
  foreach() {
    concentration[] = cs[] > 0. ? 2. : 0.;
    capacity[] = cm[];
  }
  foreach_face()
    diffusivity.x[] = fm.x[];
  boundary ({concentration});
  mgstats mg = diffusion (concentration, 0.01, D = diffusivity,
                            theta = capacity);
  double error = 0., mass = 0., wall_flux = 0.;
  foreach(reduction(max:error) reduction(+:mass) reduction(max:wall_flux)) {
    if (cs[] > 0.)
      error = max(error, fabs(concentration[] - 2.));
    mass += concentration[]*cm[]*sq(Delta);
    double flux;
    double coefficient = embed_flux(point, concentration, diffusivity, &flux);
    wall_flux = max(wall_flux, fabs(flux + coefficient*concentration[]));
  }
  double exact_mass = 2.*L0*(TEST_PIPE ? wall*wall/2. : 2.*wall);
  int bad = !isfinite(mass) || !isfinite(mg.resa) || error > 1e-9 ||
    fabs(mass - exact_mass) > 1e-9 || wall_flux > 1e-12;
  fprintf(stderr, "diffusion pipe=%d cycle=%d error=%g mass_error=%g "
          "wall_flux=%g residual=%g\n", TEST_PIPE, cycle, error,
          fabs(mass - exact_mass), wall_flux, mg.resa);
  return bad;
}

/**
A positive, nonuniform scalar must retain its physical fluid-volume
integral. An otherwise identical absorbing wall is a negative control:
its scalar integral must decrease measurably.
*/
static int check_mass (int cycle)
{
  scalar capacity[];
  face vector diffusivity[];
  double initial_mass = 0.;
  foreach(reduction(+:initial_mass)) {
    double value = cs[] > 0. ?
      2. + 0.2*cos(pi*x/2.) + 0.1*cos(pi*y/wall) : 0.;
    concentration[] = absorbing[] = value;
    initial_mass += value*cm[]*sq(Delta);
  }
  foreach_face()
    diffusivity.x[] = fm.x[];
  boundary({concentration, absorbing});
  int bad = 0;
  for (int step = 0; step < 3; step++) {
    for (scalar field in {concentration, absorbing}) {
      foreach()
        capacity[] = cm[];
      mgstats mg = diffusion(field, 0.01, D = diffusivity, theta = capacity);
      bad |= !isfinite(mg.resa) || mg.resa > 1e-9;
    }
  }
  double closed_mass = 0., absorbing_mass = 0., change = 0.;
  foreach(reduction(+:closed_mass) reduction(+:absorbing_mass)
          reduction(max:change)) {
    closed_mass += concentration[]*cm[]*sq(Delta);
    absorbing_mass += absorbing[]*cm[]*sq(Delta);
    if (cs[] > 0.) {
      double initial = 2. + 0.2*cos(pi*x/2.) + 0.1*cos(pi*y/wall);
      change = max(change, fabs(concentration[] - initial));
    }
  }
  bad |= !isfinite(closed_mass) || !isfinite(absorbing_mass) ||
    fabs(closed_mass - initial_mass) > 1e-9 ||
    initial_mass - absorbing_mass < 1e-3 || change < 1e-5;
  fprintf(stderr, "mass pipe=%d cycle=%d neumann_error=%g dirichlet_loss=%g "
          "field_change=%g\n", TEST_PIPE, cycle,
          fabs(closed_mass - initial_mass), initial_mass - absorbing_mass, change);
  return bad;
}

int main()
{
  size(4.);
  origin(-2., TEST_PIPE ? 0. : -2.);
  init_grid(16);
  TOLERANCE = 1e-11;
  run();
  return failed;
}

event init (i = 0)
{
  concentration.refine = refine_embed_linear;
  set_prolongation(concentration, refine_embed_linear);
  set_restriction(concentration, restriction_volume_average);
  absorbing.refine = refine_embed_linear;
  set_prolongation(absorbing, refine_embed_linear);
  set_restriction(absorbing, restriction_volume_average);
  confined_geometry(wall);
  for (int cycle = 0; cycle < 4; cycle++) {
    unrefine(level > 4);
    confined_geometry(wall);
    double centre = cycle % 2 ? 0.75 : -0.75;
    refine(level < 7 && fabs(x - centre) < 0.55 &&
           fabs(y - wall) < 0.3);
    confined_geometry(wall);
    failed |= check_geometry(cycle);
    failed |= check_diffusion(cycle);
    failed |= check_mass(cycle);
  }
  fprintf(stderr, "%s: embedded %s contracts\n", failed ? "FAIL" : "PASS",
          TEST_PIPE ? "pipe" : "channel");
  return 1;
}
