/**
# dropMove.c: onset of self-propulsion of a chemically active drop

A single planar drop emits a chemical species at its interface. The species
is advected and diffuses in the outer phase and lowers the interfacial
tension where it accumulates. Above a critical Péclet number the isotropic
state is unstable and the drop self-propels. The problem is integrated in
the Stokes limit with the coupled level-set and volume-of-fluid (CLSVOF)
interface method and the integral formulation of surface tension.

## Non-dimensional groups

- `Oh`: Ohnesorge number, entering only through the density $\rho = 4/Oh^2$.
- `Pe`: Péclet number; the species diffusivity is $D = 1/Pe$.
- `Ca`: the clean-interface surface tension is $1/Ca$; the local coefficient
  is $1/Ca + 4c$ with $c$ the species concentration.
- `AcNum`: constant chemical flux emitted at the interface.

## Runtime parameters

All parameters are read from a `key=value` file passed as the first
argument (default `case.params`) through `src-local/params.h`:

| Key | Default | Meaning |
|---|---|---|
| `CaseNo` | 1000 | Case identifier used by the runner |
| `Pe` | 1.6 | Péclet number |
| `MAXlevel`, `MINlevel` | 8, 0 | Quadtree refinement bounds |
| `Oh`, `Ca`, `AcNum` | 1, 0.1, 1 | Fixed non-dimensional groups |
| `L0` | 10 | Square domain size in drop radii |
| `tmax`, `tsnap` | 50, 0.1 | Observation horizon and snapshot interval |
| `threshold` | 1 | Centroid displacement, in drop radii, classified as `MOVED`; `<= 0` disables the early stop |
| `FErr`, `VelErr`, `cErr`, `KErr` | 1e-3 each | Wavelet adaptation tolerances |
| `keLimit` | 1e3 | Kinetic-energy limit classified as `FAILED` |

## Output and classification

Snapshots are written to `intermediate/snapshot-<t>` and the diagnostics to
`log.dat` with columns `i t ke dist xcm ycm volume` in the working
directory (the runner uses `simulationCases/c<CaseNo>/`). Exactly one
`STATUS {MOVED,NOT_MOVED,FAILED}` line and one `SUMMARY key=value` line are
printed on stdout when the run ends. `MOVED` is a finite-time,
finite-displacement convention: the centroid moved by more than
`threshold` before `tmax`. The centroid, drop volume and kinetic energy are
integrated with the cell volume `dv()` so that the adaptive mesh does not
bias the classifier; displacement is measured from the centroid recorded at
the first step.

## Author
Vatsal Sanjay
Email: vatsal.sanjay@comphy-lab.org
Computational Multiphase Physics (CoMPhy) Lab, Durham University
Last updated: Sep 9, 2026
*/

#include <sys/stat.h>
#include <errno.h>

#include "navier-stokes/centered.h"
#define FILTERED 1
#include "two-phase-clsvof.h"
#include "integral.h"
#include "activity.h"
#include "params.h"

/**
## Fields and boundary conditions

`cL` is the species concentration in the outer phase (`cL.inverse = true`),
`sigmaf` the surface tension coefficient and `KAPPA` the distance-function
curvature used for adaptation. All walls are no-slip with zero
concentration.
*/
scalar cL[], * stracers = {cL};
scalar sigmaf[];
scalar KAPPA[];

cL[top] = dirichlet(0.);
cL[right] = dirichlet(0.);
cL[left] = dirichlet(0.);
cL[bottom] = dirichlet(0.);

u.t[top] = dirichlet(0.);
u.t[right] = dirichlet(0.);
u.t[left] = dirichlet(0.);
u.t[bottom] = dirichlet(0.);

/**
## Runtime parameters
*/
int CaseNo = 1000;
double Pe = 1.6;
int MAXlevel = 8, MINlevel = 0;
double Oh = 1., Ca = 0.1, AcNum = 1.;
double tmax = 50., tsnap = 0.1;
double movement_threshold = 1.0;
double FErr = 1e-3, VelErr = 1e-3, cErr = 1e-3, KErr = 1e-3;
double keLimit = 1e3;

/**
## Classifier state
*/
static double x0_cm = 0., y0_cm = 0.;
static bool centroid_initialised = false;
static bool status_printed = false;
static int exit_status = 0;
static double dist_last = 0., xcm_last = 0., ycm_last = 0.;

/**
### print_status()

Prints the single `STATUS` line and the `SUMMARY` line on rank zero, once.
*/
static void print_status (const char * status, int i_now, double t_now)
{
  if (status_printed || pid() != 0)
    return;
  status_printed = true;
  fprintf (stdout, "STATUS %s\n", status);
  fprintf (stdout,
           "SUMMARY CaseNo=%d Pe=%g max_level=%d tmax=%g threshold=%g"
           " t_end=%g i_end=%d dist_end=%.8e xcm_end=%.8e ycm_end=%.8e"
           " status=%s\n",
           CaseNo, Pe, MAXlevel, tmax, movement_threshold,
           t_now, i_now, dist_last, xcm_last, ycm_last, status);
  fflush (stdout);
}

/**
### ensure_directory()

Creates `path` with `mkdir(2)`; an existing directory is accepted.
*/
static bool ensure_directory (const char * path)
{
  if (mkdir (path, 0775) != 0 && errno != EEXIST)
    return false;
  struct stat st;
  return stat (path, &st) == 0 && S_ISDIR (st.st_mode);
}

/**
## main()

Reads the parameter file, validates the values that must be positive, sets
up the domain and material properties and runs the simulation.
*/
int main (int argc, char const * argv[])
{
  params_init_from_argv (argc, argv);

  CaseNo = param_int ("CaseNo", CaseNo);
  Pe = param_double ("Pe", Pe);
  MAXlevel = param_int ("MAXlevel", MAXlevel);
  MINlevel = param_int ("MINlevel", MINlevel);
  Oh = param_double ("Oh", Oh);
  Ca = param_double ("Ca", Ca);
  AcNum = param_double ("AcNum", AcNum);
  double L0_param = param_double ("L0", 10.);
  tmax = param_double ("tmax", tmax);
  tsnap = param_double ("tsnap", tsnap);
  movement_threshold = param_double ("threshold", movement_threshold);
  FErr = param_double ("FErr", FErr);
  VelErr = param_double ("VelErr", VelErr);
  cErr = param_double ("cErr", cErr);
  KErr = param_double ("KErr", KErr);
  keLimit = param_double ("keLimit", keLimit);

  if (!(Pe > 0.) || !(Oh > 0.) || !(Ca > 0.) || !(L0_param > 0.) ||
      !(tmax > 0.) || !(tsnap > 0.) || !(keLimit > 0.) ||
      MAXlevel < 1 || MAXlevel > 20 || MINlevel < 0 || MINlevel > MAXlevel) {
    fprintf (stderr, "Invalid parameters: Pe, Oh, Ca, L0, tmax, tsnap and keLimit"
             " must be positive and 0 <= MINlevel <= MAXlevel <= 20.\n");
    return 2;
  }

  if (pid() == 0 && !ensure_directory ("intermediate")) {
    fprintf (stderr, "Cannot create output directory 'intermediate'.\n");
    return 2;
  }

  stokes = true;
  L0 = L0_param;
  origin (-0.5*L0, -0.5*L0);
  init_grid (1 << MAXlevel);

  d.sigmaf = sigmaf;

  rho1 = 4./sq(Oh); rho2 = 4./sq(Oh);
  mu1 = 1.0; mu2 = 1.0;

  cL.inverse = true;
  cL.A = AcNum;
  cL.D = 1./Pe;

  if (pid() == 0)
    fprintf (stderr, "CaseNo=%d Pe=%g MAXlevel=%d MINlevel=%d Oh=%g Ca=%g AcNum=%g"
             " L0=%g tmax=%g tsnap=%g threshold=%g\n",
             CaseNo, Pe, MAXlevel, MINlevel, Oh, Ca, AcNum, L0, tmax, tsnap,
             movement_threshold);

  run();
  return exit_status;
}

/**
## Initialisation

Signed distance to a unit circle centred at the origin; the CLSVOF header
derives the volume fraction from `d`.
*/
event init (i = 0) {
  foreach() {
    d[] = 1. - sqrt (sq(x) + sq(y));
    u.x[] = 0.0;
    u.y[] = 0.0;
    cL[] = 0.;
    sigmaf[] = 1./Ca + 4.*cL[];
  }
}

/**
## Surface tension coefficient from the local concentration
*/
event properties (i++) {
  foreach()
    sigmaf[] = 1./Ca + 4.*cL[];
}

/**
## Adaptation

Wavelet adaptation on the volume fraction, velocity, concentration and the
distance-function curvature.
*/
event adapt (i++) {
  foreach()
    KAPPA[] = distance_curvature (point, d);
  adapt_wavelet ({f, u.x, u.y, cL, KAPPA},
                 (double[]){FErr, VelErr, VelErr, cErr, KErr},
                 MAXlevel, MINlevel);
}

/**
## Snapshots
*/
event outputs (t = 0.; t += tsnap; t <= tmax) {
  char dumpFile[128];
  snprintf (dumpFile, sizeof(dumpFile), "intermediate/snapshot-%5.4f", t);
  dump (file = dumpFile);
}

/**
## Diagnostics and onset classifier

Kinetic energy, drop volume and the volume-weighted centroid are reduced
with `dv()`. The run stops with `STATUS FAILED` if the drop volume
vanishes or the kinetic energy is non-finite or exceeds `keLimit`, and with
`STATUS MOVED` once the displacement from the initial centroid exceeds
`threshold` (when `threshold > 0`).
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
    static bool log_header_written = false;
    if (!log_header_written) {
      log_header_written = true;
      fp = fopen ("log.dat", "w");
      if (fp)
        fprintf (fp, "i t ke dist xcm ycm volume\n");
      else
        fprintf (stderr, "Cannot open log.dat for writing.\n");
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

  if (!isfinite (ke) || (i > 10 && ke >= keLimit)) {
    if (pid() == 0)
      fprintf (stderr, "Kinetic energy %.8e at i=%d, t=%g: stopping.\n",
               ke, i, t);
    print_status ("FAILED", i, t);
    exit_status = 3;
    return 1;
  }

  if (movement_threshold > 0. && i > 10 && dist_last >= movement_threshold) {
    print_status ("MOVED", i, t);
    return 1;
  }
}

/**
## End of the observation window
*/
event end (t = tmax) {
  print_status ("NOT_MOVED", i, t);
}
