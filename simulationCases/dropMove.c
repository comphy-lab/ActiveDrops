/**
# dropMove.c: onset of self-propulsion of a chemically active drop

A single planar drop emits a coarse-grained chemical product at its interface.
The product is advected and diffuses in the outer phase. The implemented law
`sigma = 1/Ca + GammaSlope*cL` increases surface tension with concentration.
Above a critical Péclet number the isotropic
state is unstable and the drop self-propels. The problem uses the full
dimensionless Navier--Stokes equations with the coupled level-set and volume-of-fluid (CLSVOF)
interface method and the integral formulation of surface tension.

## Non-dimensional groups

- `Re`: outer-fluid Reynolds number; the outer density is `Re` when the outer
  viscosity, velocity and radius scales are one.
- `Ca`: clean-interface capillary number; the reference tension is $1/Ca$.
- `Pe`: product Péclet number; its diffusivity is $1/Pe$.
- `GammaSlope`: dimensionless surface-tension sensitivity; the local tension
  is $1/Ca + \mathit{GammaSlope}\,c$.
- `AcNum`: dimensionless outward normal-gradient magnitude. The PLIC
  surface-source coefficient is `AcNum/Pe`.
- `viscosityRatio`, `densityRatio`: inner-to-outer material-property ratios.

The reference-tension Ohnesorge number is derived as $Oh=\sqrt{Ca/Re}$. See
[the model contract](../src-local/active-drop-model.h) for the dimensional
scales and mobility diagnostics.

## Runtime parameters

All parameters are read from a `key=value` file passed as the first
argument (default `case.params`) through `src-local/params.h`:

| Key | Default | Meaning |
|---|---|---|
| `CaseNo` | 1000 | Case identifier used by the runner |
| `Pe` | 1.6 | Péclet number |
| `MAXlevel`, `MINlevel` | 8, 0 | Quadtree refinement bounds |
| `Re`, `Ca` | 0.01, 0.1 | Inertia and reference tension |
| `GammaSlope`, `AcNum` | 4, 1 | Surface-tension coupling and interfacial gradient magnitude |
| `viscosityRatio`, `densityRatio` | 1, 1 | Inner-to-outer property ratios |
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
the first step. Periodic coordinates are unwrapped about the preceding
centroid, so `xcm` and `ycm` may leave the principal domain.

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
#include "active-drop-model.h"

/**
## Fields and boundary conditions

`cL` is the species concentration in the outer phase (`cL.inverse = true`),
`sigmaf` the surface tension coefficient and `KAPPA` the distance-function
curvature used for adaptation. Both coordinate directions are periodic for every field, including
concentration, pressure, velocity and the CLSVOF fields.
*/
scalar cL[], * stracers = {cL};
scalar sigmaf[];
scalar KAPPA[];

/**
## Runtime parameters
*/
int CaseNo = 1000;
double Re = 0.01, Ca = 0.1, Pe = 1.6;
int MAXlevel = 8, MINlevel = 0;
double GammaSlope = 4., AcNum = 1.;
double viscosityRatio = 1., densityRatio = 1.;
double OhDerived, mobilityScaleRatio, PeMobility;
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
### periodic_coordinate()

Select the periodic image nearest the previous unwrapped centroid. The
single drop must remain compact within half a period of that reference.
This keeps centroid displacement continuous across the periodic seam.
*/
static double periodic_coordinate (double coordinate, double reference)
{
  return coordinate + L0*floor((reference - coordinate)/L0 + 0.5);
}

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
           "SUMMARY CaseNo=%d Re=%g Ca=%g Pe=%g GammaSlope=%g AcNum=%g"
           " viscosityRatio=%g densityRatio=%g Oh=%g mobilityScaleRatio=%g PeMobility=%g"
           " max_level=%d tmax=%g threshold=%g"
           " t_end=%g i_end=%d dist_end=%.8e xcm_end=%.8e ycm_end=%.8e"
           " status=%s\n",
           CaseNo, Re, Ca, Pe, GammaSlope, AcNum, viscosityRatio, densityRatio,
           OhDerived, mobilityScaleRatio, PeMobility, MAXlevel, tmax, movement_threshold,
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

  if (param_present("Oh")) {
    fprintf(stderr, "Parameter 'Oh' is retired. Supply Re and Ca; "
            "the code reports Oh=sqrt(Ca/Re).\n");
    return 2;
  }

  CaseNo = param_int ("CaseNo", CaseNo);
  Re = param_double ("Re", Re);
  Pe = param_double ("Pe", Pe);
  MAXlevel = param_int ("MAXlevel", MAXlevel);
  MINlevel = param_int ("MINlevel", MINlevel);
  Ca = param_double ("Ca", Ca);
  GammaSlope = param_double ("GammaSlope", GammaSlope);
  AcNum = param_double ("AcNum", AcNum);
  viscosityRatio = param_double ("viscosityRatio", viscosityRatio);
  densityRatio = param_double ("densityRatio", densityRatio);
  double L0_param = param_double ("L0", 10.);
  tmax = param_double ("tmax", tmax);
  tsnap = param_double ("tsnap", tsnap);
  movement_threshold = param_double ("threshold", movement_threshold);
  FErr = param_double ("FErr", FErr);
  VelErr = param_double ("VelErr", VelErr);
  cErr = param_double ("cErr", cErr);
  KErr = param_double ("KErr", KErr);
  keLimit = param_double ("keLimit", keLimit);

  if (!(Re > 0.) || !(Pe > 0.) || !(Ca > 0.) || !(L0_param > 0.) ||
      !(GammaSlope >= 0.) || !(AcNum >= 0.) ||
      !(viscosityRatio > 0.) || !(densityRatio > 0.) ||
      !(tmax > 0.) || !(tsnap > 0.) || !(keLimit > 0.) ||
      !isfinite(Re) || !isfinite(Pe) || !isfinite(Ca) ||
      !isfinite(GammaSlope) || !isfinite(AcNum) ||
      !isfinite(viscosityRatio) || !isfinite(densityRatio) ||
      MAXlevel < 1 || MAXlevel > 20 || MINlevel < 0 || MINlevel > MAXlevel) {
    fprintf (stderr, "Invalid dimensionless parameters or refinement levels. "
             "Require Re, Ca, Pe, viscosityRatio and densityRatio > 0; "
             "GammaSlope and AcNum >= 0; and 0 <= MINlevel <= MAXlevel <= 20.\n");
    return 2;
  }

  if (pid() == 0 && !ensure_directory ("intermediate")) {
    fprintf (stderr, "Cannot create output directory 'intermediate'.\n");
    return 2;
  }

  stokes = false;
  L0 = L0_param;
  origin (-0.5*L0, -0.5*L0);
  periodic (right);
  periodic (top);
  init_grid (1 << MAXlevel);

  d.sigmaf = sigmaf;

  rho1 = densityRatio*Re; rho2 = Re;
  mu1 = viscosityRatio;   mu2 = 1.;

  cL.inverse = true;
  cL.A = AcNum/Pe;
  cL.D = 1./Pe;

  OhDerived = active_drop_ohnesorge(Re, Ca);
  mobilityScaleRatio = active_drop_velocity_scale_ratio
    (AcNum, GammaSlope, viscosityRatio, false);
  PeMobility = mobilityScaleRatio*Pe;

  if (pid() == 0)
    fprintf (stderr, "CaseNo=%d Re=%g Ca=%g Pe=%g GammaSlope=%g AcNum=%g "
             "viscosityRatio=%g densityRatio=%g Oh=%g mobilityScaleRatio=%g "
             "PeMobility=%g MAXlevel=%d MINlevel=%d L0=%g tmax=%g tsnap=%g "
             "threshold=%g\n", CaseNo, Re, Ca, Pe, GammaSlope, AcNum,
             viscosityRatio, densityRatio, OhDerived, mobilityScaleRatio, PeMobility,
             MAXlevel, MINlevel, L0, tmax, tsnap, movement_threshold);

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
    sigmaf[] = 1./Ca + GammaSlope*cL[];
  }
}

/**
## Surface tension coefficient from the local concentration
*/
event properties (i++) {
  foreach()
    sigmaf[] = 1./Ca + GammaSlope*cL[];
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
    x_moment += ff*periodic_coordinate(x, xcm_last)*dv();
    y_moment += ff*periodic_coordinate(y, ycm_last)*dv();
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
