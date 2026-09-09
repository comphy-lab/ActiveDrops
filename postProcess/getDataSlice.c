#include "navier-stokes/centered.h"

char filename[512];
int nx, ny;
double xmin, ymin, xmax, ymax;
scalar * list = NULL;
scalar f[];

int main(int a, char const *arguments[]){

  L0 = 40;
  origin(-0.5*L0, -0.5*L0);

  if (a != 8) {
    fprintf (stderr, "usage: getDataSlice <snapshot> xmin xmax ymin ymax nx ny\n");
    return 1;
  }
  if (snprintf (filename, sizeof(filename), "%s", arguments[1]) >= (int) sizeof(filename)) {
    fprintf (stderr, "error: snapshot path longer than %zu characters\n", sizeof(filename) - 1);
    return 1;
  }
  xmin = atof(arguments[2]);
  xmax = atof(arguments[3]);
  ymin = atof(arguments[4]);
  ymax = atof(arguments[5]);
  nx = atoi(arguments[6]);
  ny = atoi(arguments[7]);
  if (nx <= 0 || ny <= 0 || !(xmax > xmin) || !(ymax > ymin)) {
    fprintf (stderr, "error: need xmax > xmin, ymax > ymin and positive nx, ny\n");
    return 1;
  }

  list = list_add (list, f);
  list = list_add (list, u.x);
  list = list_add (list, u.y);

  if (!restore (file = filename)) {
    fprintf (stderr, "error: could not restore %s\n", filename);
    return 1;
  }
  boundary((scalar *){f, u.x, u.y});

  FILE * fp = ferr;
  nx++;
  ny++;
  double Deltax = 0.999999*(xmax-xmin)/(nx - 1);
  double Deltay = 0.999999*(ymax-ymin)/(ny - 1);
  int len = list_len(list);
  double ** field = (double **) matrix_new (nx, ny, len*sizeof(double));
  for (int i = 0; i < nx; i++) {
    double x = Deltax*i + xmin;
    for (int j = 0; j < ny; j++) {
      double y = Deltay*j + ymin;
      int k = 0;
      for (scalar s in list){
        field[i][len*j + k++] = interpolate (s, x, y);
      }
    }
  }

  for (int i = 0; i < nx; i++) {
    double x = Deltax*i + xmin;
    for (int j = 0; j < ny; j++) {
      double y = Deltay*j + ymin;
      fprintf (fp, "%g %g", x, y);
      int k = 0;
      for (scalar s in list){
        fprintf (fp, " %g", field[i][len*j + k++]);
      }
      fputc ('\n', fp);
    }
  }
  fflush (fp);
  matrix_free (field);

}
