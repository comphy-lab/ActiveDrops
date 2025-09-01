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

  foreach() {
    sumv1 += clamp(f[], 0., 1.)*x;
    sumv2 += clamp(f[], 0., 1.)*y;

    sumf += clamp(f[], 0., 1.);    
  }

  xcm1 = sumv1/sumf;
  ycm2 = sumv2/sumf;

  boundary((scalar *){f, u.x, u.y});

  FILE * fp = ferr;
  fprintf(fp, "%f %f %f\n", interpolate(u.x,xcm1, ycm2), interpolate(u.y,xcm1, ycm2), t);

  fflush (fp);
  fclose (fp);
}