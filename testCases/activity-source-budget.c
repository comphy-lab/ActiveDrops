/**
# Geometric activity-source budget

Check the production PLIC source on uniform grids against the analytical
circumference of a planar circle and surface area of a sphere. The
axisymmetric area omits the common azimuthal factor, consistently with the
solver metric. Pure drop cells and full solids must receive no source.

An exactly face-aligned pure-phase interface provides a separate fallback
check because it contains no mixed PLIC cell. These are bounded software and
geometric-convergence checks, not verification of the coupled transport.
*/

#include "grid/multigrid.h"
#if TEST_PIPE
# include "axi.h"
#endif
#include "run.h"
#include "vof.h"

scalar f[], * interfaces = {f};
scalar p[];
face vector uf[];
#include "activity.h"

scalar concentration[], * stracers = {concentration};
static int failed = 0;
static double previous_error = HUGE;

/**
Apply one source-only update and return the produced mass divided by
`dt*concentration.A`. The same pass checks that the production source has no
support in pure inner-phase cells.
*/
static double source_area (int * bad_support)
{
  foreach() concentration[] = 0.;
  foreach_face() uf.x[] = 0.;
  boundary({f, concentration, uf});

  dt = 0.025;
  concentration.inverse = true;
  concentration.A = 0.15;
  concentration.D = 0.;
  event("vof");
  event("tracer_diffusion");

  double mass = 0.;
  int inner_source = 0;
  foreach(reduction(+:mass) reduction(+:inner_source)) {
    mass += concentration[]*dv();
    if (f[] >= 1. - 1e-6 && ActivityFlux[] != 0.)
      inner_source++;
  }
  *bad_support = inner_source;
  return mass/(dt*concentration.A);
}

int main()
{
  size(4.);
  origin(-2., TEST_PIPE ? 0. : -2.);
  for (N = 32; N <= 128; N *= 2)
    run();
  return failed;
}

event init (i = 0)
{
  fraction(f, 1. - sqrt(sq(x) + sq(y)));
  int bad_support = 0;
  double measured = source_area(&bad_support);
  const double exact = TEST_PIPE ? 2. : 2.*pi;
  double error = fabs(measured - exact);
  int bad = !isfinite(measured) || bad_support || error > 0.2 ||
    (N > 32 && error > 1.05*previous_error);
  if (N == 128)
    bad |= error > 0.03;
  previous_error = error;
  failed |= bad;
  fprintf(stderr, "geometric-source pipe=%d N=%d area=%g exact=%g error=%g "
          "pure_inner=%d status=%s\n", TEST_PIPE, N, measured, exact, error,
          bad_support, bad ? "FAIL" : "PASS");

  foreach() f[] = x < 0. ? 1. : 0.;
  boundary({f});
  bad_support = 0;
  measured = source_area(&bad_support);
  const double face_exact = TEST_PIPE ? 8. : 4.;
  error = fabs(measured - face_exact);
  bad = !isfinite(measured) || bad_support || error > 1e-12;
  failed |= bad;
  fprintf(stderr, "face-source pipe=%d N=%d area=%g exact=%g error=%g "
          "pure_inner=%d status=%s\n", TEST_PIPE, N, measured, face_exact,
          error, bad_support, bad ? "FAIL" : "PASS");
}
