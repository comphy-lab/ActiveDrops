/**
# Activity tracer phase assignment

Exercise the production activity event with both parent concentration
phase flags. The temporary fields must represent the two complementary
VOF phases. Uniform concentration must survive a nonzero translation of
the interface with activity and diffusion disabled.
*/

#include "grid/quadtree.h"
#include "run.h"
#include "vof.h"

scalar f[], * interfaces = {f};
scalar p[];
face vector uf[];
#include "activity.h"

scalar concentration[], * stracers = {concentration};
static int failed = 0;

int main()
{
  init_grid(32);
  periodic(right);
  periodic(top);
  run();
  return failed;
}

event init (i = 0)
{
  const double uniform = 2.75;
  dt = 0.2/32.;
  concentration.A = 0.;
  concentration.D = 0.;
  scalar before[];
  for (int inverse = 0; inverse < 2; inverse++) {
    concentration.inverse = inverse;
    fraction(f, min(x - 0.213, 0.713 - x));
    foreach() {
      concentration[] = uniform;
      before[] = f[];
    }
    foreach_face()
      uf.x[] = 0.;
    boundary({concentration, f, uf});

    /* Invoke the production allocation event and a zero-velocity VOF sweep. */
    event("vof");
    scalar first = concentration.phi1, second = concentration.phi2;
    int flag_error = first.inverse || !second.inverse;
    double split_error = 0.;
    foreach(reduction(max:split_error)) {
      split_error = max(split_error, fabs(first[] - uniform*f[]));
      split_error = max(split_error, fabs(second[] - uniform*(1. - f[])));
    }

    foreach_face(x)
      uf.x[] = 1.;
    boundary({uf});
    vof_advection(interfaces, inverse);
    double sum_error = 0., displacement = 0.;
    foreach(reduction(max:sum_error) reduction(max:displacement)) {
      sum_error = max(sum_error, fabs(first[] + second[] - uniform));
      displacement = max(displacement, fabs(f[] - before[]));
    }
    event("tracer_diffusion");
    double reconstructed_error = 0.;
    foreach(reduction(max:reconstructed_error))
      reconstructed_error = max(reconstructed_error,
                                fabs(concentration[] - uniform));
    int bad = flag_error || split_error > 1e-12 || sum_error > 1e-12 ||
      reconstructed_error > 1e-12 || displacement < 0.1 || f.tracers != NULL;
    failed |= bad;
    fprintf(stderr, "activity inverse=%d flag_error=%d split_error=%g "
            "sum_error=%g reconstructed_error=%g displacement=%g\n",
            inverse, flag_error, split_error, sum_error,
            reconstructed_error, displacement);
  }
  fprintf(stderr, "%s: activity phase contracts\n", failed ? "FAIL" : "PASS");
  return 1;
}
