/* Title: Getting Facets
# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids
*/
#include "utils.h"
#include "output.h"
#include "fractions.h"

scalar f[];
char filename[512];
int main(int a, char const *arguments[])
{
  if (a != 2) {
    fprintf (stderr, "usage: getFacets <snapshot>\n");
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
  FILE * fp = ferr;
  output_facets(f,fp);
  fflush (fp);
  return 0;
}

