#include "navier-stokes/centered.h"
#include "fractions.h"

char filename[512];
double xcm1 , ycm2;

scalar f[];


int main(int a, char const *arguments[])
{
  if (a != 2) {
    fprintf (stderr, "usage: getVelocity_v2 <snapshot>\n");
    return 1;
  }
  if (snprintf (filename, sizeof(filename), "%s", arguments[1]) >= (int) sizeof(filename)) {
    fprintf (stderr, "error: snapshot path longer than %zu characters\n", sizeof(filename) - 1);
    return 1;
  }
  if (!restore (file = filename)) {
    fprintf (stderr, "error: could not restore %s\n", filename);
    return 1;
  }
  f.prolongation = fraction_refine;
  boundary((scalar *){f, u.x, u.y});

  double sumv1 = 0.;
  double sumv2 = 0.;
  double sumf = 0.;

  /* Volume-weighted centroid, consistent with getCM.c and the driver. */
  foreach (reduction(+:sumv1) reduction(+:sumv2) reduction(+:sumf)) {
    double ff = clamp(f[], 0., 1.);
    sumv1 += ff*x*dv();
    sumv2 += ff*y*dv();
    sumf += ff*dv();
  }

  if (!(sumf > 0.)) {
    fprintf (ferr, "0 0 %f\n", t);
    return 1;
  }
  xcm1 = sumv1/sumf;
  ycm2 = sumv2/sumf;

  boundary((scalar *){f, u.x, u.y});

  FILE * fp = ferr;
  fprintf(fp, "%f %f %f\n", interpolate(u.x,xcm1, ycm2), interpolate(u.y,xcm1, ycm2), t);

  fflush (fp);
  return 0;
}
