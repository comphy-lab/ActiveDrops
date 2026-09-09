/**
 * Centroid diagnostic check on a deliberately asymmetric adaptive mesh
 *
 * Two circles are initialised as volume fractions on a quadtree whose
 * right half (x > 0) is refined two levels deeper than the left half.
 * For each circle the centroid is computed twice: with cell-count
 * weights (the original dropMove diagnostic, sum f*x / sum f) and with
 * volume weights (sum f*x*dv() / sum f*dv(), as used by the fixed
 * driver). The volume-weighted centroid must recover the imposed centre
 * to within a fraction of the coarse cell size regardless of the
 * refinement; the cell-count centroid is biased towards the refined half.
 *
 * Compile and run from the testCases directory (or use testCases/run-tests.sh):
 *   qcc -O2 -Wall -disable-dimensions centroid-check.c -o centroid-check -lm
 *   ./centroid-check
 * Exit status 0 means both volume-weighted checks passed.
 */

#include "grid/quadtree.h"
#include "fractions.h"
#include "utils.h"

#define BASE_LEVEL 6
#define FINE_LEVEL 8

scalar f[];

static void centroids (double * xc_count, double * yc_count,
                       double * xc_vol, double * yc_vol, double * vol)
{
  double sf = 0., sx = 0., sy = 0., vf = 0., vx = 0., vy = 0.;
  foreach (reduction(+:sf) reduction(+:sx) reduction(+:sy)
           reduction(+:vf) reduction(+:vx) reduction(+:vy)) {
    double ff = clamp (f[], 0., 1.);
    sf += ff; sx += ff*x; sy += ff*y;
    vf += ff*dv(); vx += ff*x*dv(); vy += ff*y*dv();
  }
  *xc_count = sx/sf; *yc_count = sy/sf;
  *xc_vol = vx/vf; *yc_vol = vy/vf;
  *vol = vf;
}

static int check_circle (double cx, double cy, double radius, const char * label)
{
  refine (x > 0. && level < FINE_LEVEL);
  vertex scalar phi[];
  foreach_vertex()
    phi[] = sq(radius) - sq(x - cx) - sq(y - cy);
  fractions (phi, f);

  double xc_count, yc_count, xc_vol, yc_vol, vol;
  centroids (&xc_count, &yc_count, &xc_vol, &yc_vol, &vol);
  double coarse = L0/(1 << BASE_LEVEL);
  double err_vol = sqrt (sq(xc_vol - cx) + sq(yc_vol - cy));
  double err_count = sqrt (sq(xc_count - cx) + sq(yc_count - cy));
  double err_area = fabs (vol - pi*sq(radius))/(pi*sq(radius));
  int ok = err_vol < 0.05*coarse && err_area < 1e-2;
  printf ("%s: centre=(%g,%g) r=%g\n", label, cx, cy, radius);
  printf ("  cell-count centroid   = (%.6e, %.6e)  error %.3e (%.2f coarse cells)\n",
          xc_count, yc_count, err_count, err_count/coarse);
  printf ("  volume-weighted centroid = (%.6e, %.6e)  error %.3e (%.4f coarse cells)\n",
          xc_vol, yc_vol, err_vol, err_vol/coarse);
  printf ("  volume = %.8e (relative error %.3e)  -> %s\n", vol, err_area,
          ok ? "PASS" : "FAIL");
  unrefine (level > BASE_LEVEL);
  return ok;
}

int main()
{
  L0 = 10.;
  origin (-0.5*L0, -0.5*L0);
  init_grid (1 << BASE_LEVEL);

  int ok = 1;
  ok &= check_circle (0., 0., 1., "stationary drop at the origin");
  ok &= check_circle (0.3, -0.2, 1., "translated drop");
  printf ("%s\n", ok ? "ALL PASS" : "SOME FAIL");
  return ok ? 0 : 1;
}
