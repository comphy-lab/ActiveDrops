#include "navier-stokes/centered.h"
#include "fractions.h"

char filename[80];
double xcm1 , ycm2;

scalar f[];


int main(int a, char const *arguments[])
{
  sprintf (filename, "%s", arguments[1]);

  restore (file = filename);
  f.prolongation = fraction_refine;
  boundary((scalar *){f, u.x, u.y});

  double sumv1 = 0.;
  double sumv2 = 0.;
  double sumf = 0.;

  /* Volume-weighted centroid: on an adaptive mesh the cell count is not a
     measure of area, so each cell contributes f*dv(). */
  foreach (reduction(+:sumv1) reduction(+:sumv2) reduction(+:sumf)) {
    double ff = clamp(f[], 0., 1.);
    sumv1 += ff*x*dv();
    sumv2 += ff*y*dv();
    sumf += ff*dv();
  }

  xcm1 = sumf > 0. ? sumv1/sumf : 0.;
  ycm2 = sumf > 0. ? sumv2/sumf : 0.;

  boundary((scalar *){f, u.x, u.y});

  FILE * fp = ferr;
  fprintf(fp, "%f %f %f\n", xcm1, ycm2, t);

  fflush (fp);
  fclose (fp);
}