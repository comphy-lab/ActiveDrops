/* Title: getting Data from simulation snapshot
# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids
*/

#include "utils.h"
#include "output.h"

vector u[];

char filename[512];
int nx, ny, len;
double xmin, ymin, xmax, ymax, Deltax, Deltay;
scalar * list = NULL;
scalar cL[], D2c[], f[], vel[];

int main(int a, char const *arguments[])
{
  if (a != 7) {
    fprintf (stderr, "usage: getData <snapshot> xmin ymin xmax ymax ny\n");
    return 1;
  }
  if (snprintf (filename, sizeof(filename), "%s", arguments[1]) >= (int) sizeof(filename)) {
    fprintf (stderr, "error: snapshot path longer than %zu characters\n", sizeof(filename) - 1);
    return 1;
  }
  xmin = atof(arguments[2]); ymin = atof(arguments[3]);
  xmax = atof(arguments[4]); ymax = atof(arguments[5]);
  ny = atoi(arguments[6]);
  if (ny <= 0 || !(xmax > xmin) || !(ymax > ymin)) {
    fprintf (stderr, "error: need xmax > xmin, ymax > ymin and positive ny\n");
    return 1;
  }

  list = list_add (list, cL);
  list = list_add (list, D2c);
  list = list_add (list, vel);

  /*
  Actual run and codes!
  */
  if (!restore (file = filename)) {
    fprintf (stderr, "error: could not restore %s\n", filename);
    return 1;
  }

  foreach(){
    double ff = clamp(f[], 0., 1.);

    cL[] *= (1.-ff);
    vel[] = sqrt(sq(u.x[]) + sq(u.y[]));

    double D11 = (u.y[0,1] - u.y[0,-1])/(2*Delta);
    double D22 = (u.x[1,0] - u.x[-1,0])/(2*Delta);
    double D12 = 0.5*( (u.y[1,0] - u.y[-1,0] + u.x[0,1] - u.x[0,-1])/(2*Delta) );
    double D2 = sqrt(sq(D11)+sq(D22)+2.0*sq(D12));
    D2c[] = D2/sqrt(2.0);

    if (D2c[] > 0.){
      D2c[] = log(D2c[])/log(10);
    } else {
      D2c[] = -10;
    }
    
  }

  FILE * fp = ferr;
  Deltay = (double)((ymax-ymin)/(ny));
  // fprintf(ferr, "%g\n", Deltay);
  nx = (int)((xmax - xmin)/Deltay);
  // fprintf(ferr, "%d\n", nx);
  Deltax = (double)((xmax-xmin)/(nx));
  // fprintf(ferr, "%g\n", Deltax);
  len = list_len(list);
  // fprintf(ferr, "%d\n", len);
  double ** field = (double **) matrix_new (nx, ny+1, len*sizeof(double));
  for (int i = 0; i < nx; i++) {
    double x = Deltax*(i+1./2) + xmin;
    for (int j = 0; j < ny; j++) {
      double y = Deltay*(j+1./2) + ymin;
      int k = 0;
      for (scalar s in list){
        field[i][len*j + k++] = interpolate (s, x, y);
      }
    }
  }

  for (int i = 0; i < nx; i++) {
    double x = Deltax*(i+1./2) + xmin;
    for (int j = 0; j < ny; j++) {
      double y = Deltay*(j+1./2) + ymin;
      fprintf (fp, "%g %g", x, y);
      int k = 0;
      for (scalar s in list){
        fprintf (fp, " %g", field[i][len*j + k++]);
      }
      fputc ('\n', fp);
    }
  }
  fflush (fp);
  return 0;
  matrix_free (field);
}
