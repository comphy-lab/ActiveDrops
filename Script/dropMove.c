/**
 * Onset of self-propulsion of a chemically active drop: Péclet-number scan driver
 *
 * A single planar drop with an interfacial chemical source and a
 * concentration-dependent surface tension is integrated in the Stokes
 * limit using the CLSVOF (coupled level-set and volume-of-fluid) method
 * and the integral formulation of surface tension.
 *
 * Non-dimensional groups (see README.md):
 * - Oh (Ohnesorge number): sets the density through rho = 4/Oh^2
 * - Pe (Péclet number): ratio of advection to diffusion of the chemical
 * - Ca: inverse surface tension coefficient (lower Ca, more circular drop)
 * - AcNum: constant chemical flux emitted at the drop interface
 *
 * Usage:
 *   ./dropMove <Pe> [tmax=50] [tsnap=0.1] [max_level=9] [threshold=1] [out=intermediate]
 *
 * The run terminates when the drop centroid has moved by more than
 * `threshold` (in units of the drop radius) from its initial centroid, or
 * when t reaches tmax. Exactly one classification line is printed on
 * stdout, `STATUS MOVED`, `STATUS NOT_MOVED` or `STATUS FAILED`, followed
 * by a single `SUMMARY` line recording the run parameters, the final
 * time and the final displacement. The displacement history is written to
 * `<out>/log.dat` so that growth rates near onset can be extracted.
 *
 * The centroid is computed with volume (area) weights, `f*dv()`, so that
 * the adaptive mesh does not bias the classifier.
 */

#define MIN_LEVEL 0

#define VelErr 1e-3
#define FErr 1e-3
#define cErr 1e-3
#define KErr 1e-3

#include <sys/stat.h>
#include <errno.h>

#include "navier-stokes/centered.h"
#define FILTERED 1
#include "two-phase-clsvof.h"
#include "integral.h"
#include "../src-local/activity.h"

/**
 * Global variables and boundary conditions
 * cL: chemical concentration field (in the outer phase, cL.inverse = true)
 * sigmaf: surface tension coefficient field
 */
scalar cL[],  *stracers = {cL};
#define c0 0.0

cL[top] = dirichlet(0.);
cL[right] = dirichlet(0.);
cL[left] = dirichlet(0.);
cL[bottom] = dirichlet(0.);

u.t[top] = dirichlet(0.);
u.t[right] = dirichlet(0.);
u.t[left] = dirichlet(0.);
u.t[bottom] = dirichlet(0.);

scalar sigmaf[];
scalar KAPPA[];

/**
 * Fixed non-dimensional parameters
 */
#define Oh 1e0
#define Ca 0.1
#define AcNum 1e0

/**
 * Runtime parameters
 */
double Pe = 0.;
double tmax = 50.;
double tsnap = 1e-1;
int max_level = 9;
double movement_threshold = 1.0;
char output_dir[256] = "intermediate";

/**
 * Runtime state used by the onset classifier
 */
static double x0_cm = 0., y0_cm = 0.;   // initial centroid
static bool centroid_initialised = false;
static bool status_printed = false;
static int exit_status = 0;
static double dist_last = 0., xcm_last = 0., ycm_last = 0.;

static void print_status (const char * status, int i_now, double t_now)
{
  if (status_printed || pid() != 0)
    return;
  status_printed = true;
  fprintf (stdout, "STATUS %s\n", status);
  fprintf (stdout,
           "SUMMARY Pe=%g max_level=%d tmax=%g threshold=%g"
           " t_end=%g i_end=%d dist_end=%.8e xcm_end=%.8e ycm_end=%.8e"
           " status=%s\n",
           Pe, max_level, tmax, movement_threshold,
           t_now, i_now, dist_last, xcm_last, ycm_last, status);
  fflush (stdout);
}

static bool parse_double (const char * value, double * out)
{
  char * end = NULL;
  double v = strtod (value, &end);
  if (end == value || *end != '\0' || !isfinite (v))
    return false;
  *out = v;
  return true;
}

/**
 * Create `path` and any missing parent components with mkdir(2); an
 * existing directory is accepted.
 */
static bool ensure_directory (const char * path)
{
  char buffer[sizeof(output_dir)];
  size_t n = strlen (path);
  if (n == 0 || n >= sizeof(buffer))
    return false;
  memcpy (buffer, path, n + 1);
  for (char * p = buffer + 1; *p; p++)
    if (*p == '/') {
      *p = '\0';
      if (mkdir (buffer, 0775) != 0 && errno != EEXIST)
        return false;
      *p = '/';
    }
  if (mkdir (buffer, 0775) != 0 && errno != EEXIST)
    return false;
  struct stat st;
  return stat (buffer, &st) == 0 && S_ISDIR (st.st_mode);
}

static void usage (const char * prog)
{
  fprintf (stderr,
           "usage: %s <Pe> [tmax=50] [tsnap=0.1] [max_level=9]"
           " [threshold=1] [out=intermediate]\n", prog);
}

static bool parse_arguments (int argc, char const * argv[])
{
  if (argc < 2 || !parse_double (argv[1], &Pe) || Pe <= 0.) {
    usage (argv[0]);
    return false;
  }
  for (int k = 2; k < argc; k++) {
    char key[64], value[256];
    const char * eq = strchr (argv[k], '=');
    if (!eq || eq == argv[k] || (size_t)(eq - argv[k]) >= sizeof(key) ||
        strlen (eq + 1) >= sizeof(value)) {
      fprintf (stderr, "Cannot parse argument '%s'.\n", argv[k]);
      usage (argv[0]);
      return false;
    }
    memcpy (key, argv[k], eq - argv[k]);
    key[eq - argv[k]] = '\0';
    strcpy (value, eq + 1);
    double v;
    if (!strcmp (key, "tmax") && parse_double (value, &v) && v > 0.)
      tmax = v;
    else if (!strcmp (key, "tsnap") && parse_double (value, &v) && v > 0.)
      tsnap = v;
    else if (!strcmp (key, "max_level") && parse_double (value, &v) &&
             v >= 1. && v <= 20. && v == (int) v)
      max_level = (int) v;
    else if (!strcmp (key, "threshold") && parse_double (value, &v) && v > 0.)
      movement_threshold = v;
    else if (!strcmp (key, "out"))
      strcpy (output_dir, value);
    else {
      fprintf (stderr, "Invalid argument '%s'.\n", argv[k]);
      usage (argv[0]);
      return false;
    }
  }
  return true;
}

/**
 * Main function: sets up and runs the simulation
 */
int main (int argc, char const * argv[])
{
  if (!parse_arguments (argc, argv))
    return 2;

  stokes = true;
  L0 = 10.0;
  origin (-0.5*L0, -0.5*L0);
  N = 1 << max_level;
  init_grid (N);

  d.sigmaf = sigmaf;

  rho1 = 4./sq(Oh); rho2 = 4./sq(Oh);
  mu1 = 1.0; mu2 = 1.0;

  cL.inverse = true;
  cL.A = AcNum;
  cL.D = 1./Pe;

  if (pid() == 0 && !ensure_directory (output_dir)) {
    fprintf (stderr, "Cannot create output directory '%s'.\n", output_dir);
    return 2;
  }

  run();
  return exit_status;
}

/**
 * Initialisation: signed distance to a unit circle centred at the origin;
 * the CLSVOF header derives the volume fraction from d.
 */
event init (i = 0) {
  foreach() {
    d[] = 1. - sqrt (sq(x) + sq(y));
    u.x[] = 0.0;
    u.y[] = 0.0;
    cL[] = c0;
    sigmaf[] = 1./Ca + 4.*cL[];
  }
}

/**
 * Surface tension coefficient from the local concentration
 */
event properties (i++) {
  foreach()
    sigmaf[] = 1./Ca + 4.*cL[];
}

event adapt (i++) {
  foreach()
    KAPPA[] = distance_curvature (point, d);
  adapt_wavelet ({f, u.x, u.y, cL, KAPPA},
                 (double[]){FErr, VelErr, VelErr, cErr, KErr},
                 max_level, MIN_LEVEL);
}

/**
 * Snapshots at regular intervals
 */
event outputs (t = 0.; t += tsnap; t <= tmax) {
  char dumpFile[320];
  snprintf (dumpFile, sizeof(dumpFile), "%s/snapshot-%5.4f", output_dir, t);
  dump (file = dumpFile);
}

/**
 * Diagnostics and onset classifier
 *
 * The kinetic energy and the drop centroid are both integrated with the
 * cell volume dv(), so refinement does not change their value. The
 * displacement is measured from the centroid recorded at i = 0.
 */
event logWriting (i++) {
  double ke = 0., drop_volume = 0., x_moment = 0., y_moment = 0.;
  foreach (reduction(+:ke) reduction(+:drop_volume)
           reduction(+:x_moment) reduction(+:y_moment)) {
    double ff = clamp (f[], 0., 1.);
    ke += 0.5*rho(ff)*(sq(u.x[]) + sq(u.y[]))*dv();
    drop_volume += ff*dv();
    x_moment += ff*x*dv();
    y_moment += ff*y*dv();
  }

  if (!(drop_volume > 0.)) {
    if (pid() == 0)
      fprintf (stderr, "Drop volume vanished at i=%d, t=%g.\n", i, t);
    print_status ("FAILED", i, t);
    exit_status = 3;
    return 1;
  }

  xcm_last = x_moment/drop_volume;
  ycm_last = y_moment/drop_volume;
  if (!centroid_initialised) {
    x0_cm = xcm_last, y0_cm = ycm_last;
    centroid_initialised = true;
  }
  dist_last = sqrt (sq(xcm_last - x0_cm) + sq(ycm_last - y0_cm));

  if (pid() == 0) {
    static FILE * fp = NULL;
    if (!fp) {
      char log_path[320];
      snprintf (log_path, sizeof(log_path), "%s/log.dat", output_dir);
      fp = fopen (log_path, "w");
      if (fp)
        fprintf (fp, "i t ke dist xcm ycm volume\n");
      fprintf (ferr, "i t ke dist xcm ycm volume\n");
    }
    fprintf (ferr, "%d %g %.8e %.8e %.8e %.8e %.8e\n",
             i, t, ke, dist_last, xcm_last, ycm_last, drop_volume);
    if (fp) {
      fprintf (fp, "%d %g %.8e %.8e %.8e %.8e %.8e\n",
               i, t, ke, dist_last, xcm_last, ycm_last, drop_volume);
      fflush (fp);
    }
  }

  if (!isfinite (ke) || (i > 10 && ke >= 1e3)) {
    if (pid() == 0)
      fprintf (stderr, "Kinetic energy %.8e at i=%d, t=%g: stopping.\n",
               ke, i, t);
    print_status ("FAILED", i, t);
    exit_status = 3;
    return 1;
  }

  if (i > 10 && dist_last >= movement_threshold) {
    print_status ("MOVED", i, t);
    return 1;
  }
}

event end (t = tmax) {
  print_status ("NOT_MOVED", i, t);
}
