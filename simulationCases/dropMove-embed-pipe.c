/**
# dropMove-embed-pipe.c: chemically active drop in a pipe

A unit sphere starts on the symmetry axis `y=0`. The coordinate `x` is
axial and `y` is radial; axisymmetry excludes transverse migration.
The straight embedded pipe wall is at `y=wallHalfWidth`.
The left/right boundaries are periodic for every field. The embedded walls
are no-slip, non-wetting to the drop and impermeable to chemical species.

The Navier--Stokes/CLSVOF/integral surface-tension model is the same as
`dropMove.c`: $D=1/Pe$, $\sigma=1/Ca+\mathit{GammaSlope}\,c_L$ and geometric
PLIC source $(AcNum/Pe)\delta_{\Gamma,h}$.
Positive `AcNum` produces species, not consumed fuel. Activity adds no
solid-wall chemical flux.

## Runtime parameters

All parameters use the `key=value` file passed as `argv[1]` through
`src-local/params.h`. Lengths are in initial drop radii.

| Key | Default | Meaning |
|---|---|---|
| `CaseNo` | 1001 | Case identifier used by the runner |
| `Pe` | 1.6 | Péclet number |
| `MAXlevel`, `MINlevel` | 8, 0 | Quadtree refinement bounds |
| `Re`, `Ca` | 0.01, 0.1 | Inertia and reference tension |
| `GammaSlope`, `AcNum` | 4, 1 | Surface-tension coupling and interfacial gradient magnitude |
| `viscosityRatio`, `densityRatio` | 1, 1 | Inner-to-outer property ratios |
| `L0` | 10 | Square domain size in drop radii |
| `tmax`, `tsnap` | 50, 0.1 | Observation horizon and snapshot interval |
| `threshold` | 0 | Centroid displacement, in drop radii, classified as `MOVED`; `<= 0` disables the early stop |
| `FErr`, `VelErr`, `cErr`, `KErr` | 1e-3 each | Wavelet adaptation tolerances |
| `keLimit` | 1e3 | Kinetic-energy limit classified as `FAILED` |
| `wallHalfWidth` | 2.53 | Pipe radius |
| `dropOffset` | 0 | Must be zero to preserve axisymmetry |
| `initialDipole` | 1e-3 | Nonnegative periodic axial species seed; zero preserves fore-aft symmetry |

## Output and classification

Snapshots, `log.dat`, `STATUS` and `SUMMARY` follow `dropMove.c`.
The periodic axial centroid is unwrapped about the previous centroid,
assuming a compact single drop within half a period of that reference.
It may leave the principal domain after crossing a periodic seam.
In the pipe, volume and energy omit the common factor $2\pi$. `ycm` is
the mean radius, not a transverse centre of mass; `dist` is axial displacement.
A liquid-containing cell entering a three-finest-cell wall band stops
with `FAILED`; this case does not supply a contact-line evolution law.

## Author
Vatsal Sanjay (vatsal.sanjay@comphy-lab.org)
Computational Multiphase Physics (CoMPhy) Lab, Durham University
*/


#include <sys/stat.h>
#include <errno.h>

#include "embed.h"
#include "axi.h"
#include "navier-stokes/centered.h"
#define FILTERED 1
#include "two-phase-clsvof.h"
#include "integral.h"
#include "activity.h"
#include "params.h"
#include "embed-channel-geometry.h"
#include "active-drop-model.h"

/**
## Fields and boundary conditions

`cL` is the species concentration in the outer phase (`cL.inverse = true`),
`sigmaf` the surface tension coefficient and `KAPPA` the distance-function
curvature used for adaptation. The embedded walls are no-slip and
non-wetting to the dispersed phase: `f=0` at the wall and in full solids,
with a negative distance field in the solid. `d[embed]` supplies numerical
extrapolation, not a prescribed contact angle. Species flux is zero
independently of the phase-wall condition.
Left/right boundaries are periodic. The bottom is the symmetry axis.
*/
scalar cL[], * stracers = {cL};
scalar sigmaf[];
scalar KAPPA[];

cL[embed] = neumann(0.);
cL[top] = neumann(0.);
cL[bottom] = neumann(0.);

u.n[top] = dirichlet(0.);
u.t[top] = dirichlet(0.);
u.n[embed] = dirichlet(0.);
u.t[embed] = dirichlet(0.);
d[embed] = neumann(0.);
f[embed] = dirichlet(0.);
f[top] = dirichlet(0.);

/**
## Runtime parameters
*/
int CaseNo = 1001;
double Re = 0.01, Ca = 0.1, Pe = 1.6;
int MAXlevel = 8, MINlevel = 0;
double GammaSlope = 4., AcNum = 1.;
double viscosityRatio = 1., densityRatio = 1.;
double OhDerived, mobilityScaleRatio, PeMobility;
double tmax = 50., tsnap = 0.1;
double movement_threshold = 0.;
double wallHalfWidth = 2.53, initialDipole = 1e-3;
double dropOffset = 0.;
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
  wallHalfWidth = param_double ("wallHalfWidth", wallHalfWidth);
  dropOffset = param_double ("dropOffset", dropOffset);
  initialDipole = param_double ("initialDipole", initialDipole);
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

  double dx = L0_param/(1 << MAXlevel);
  // A face-aligned wall has no cut-cell fragment for embedded no-slip flux.
  double wall_index = wallHalfWidth/dx;
  if (isfinite(wall_index) && fabs(wall_index - round(wall_index)) < 1e-8) {
    fprintf(stderr, "wallHalfWidth places a wall on a finest-grid face. "
            "Choose a nonaligned value (default 2.53 for L0=10).\n");
    return 2;
  }
  if (!isfinite(wallHalfWidth) || !isfinite(dropOffset) ||
      !isfinite(initialDipole) || initialDipole < 0. ||
      !isfinite(L0_param) ||
      !(wallHalfWidth > 1. + fabs(dropOffset) + 4.*dx) ||
      !(L0_param/2. > 1. + 4.*dx) ||
      dropOffset != 0. || !(wallHalfWidth < L0_param - 2.*dx) ||
      !(FErr > 0.) || !(VelErr > 0.) || !(cErr > 0.) || !(KErr > 0.) ||
      !isfinite(Pe) || !isfinite(Ca) ||
      !isfinite(tmax) || !isfinite(tsnap) || !isfinite(keLimit) ||
      !isfinite(movement_threshold) || !isfinite(FErr) ||
      !isfinite(VelErr) || !isfinite(cErr) || !isfinite(KErr)) {
    fprintf(stderr, "Invalid embedded geometry, dipole or adaptation tolerance. "
            "Keep the unit drop at least four finest cells from walls and the initial periodic-image midpoint; "
            "pipe dropOffset must be zero.\n");
    return 2;
  }

  if (pid() == 0 && !ensure_directory ("intermediate")) {
    fprintf (stderr, "Cannot create output directory 'intermediate'.\n");
    return 2;
  }

  stokes = false;
  L0 = L0_param;
  origin (-0.5*L0, 0.);
  periodic (right);
  init_grid (1 << MAXlevel);

  d.sigmaf = sigmaf;

  rho1 = densityRatio*Re; rho2 = Re;
  mu1 = viscosityRatio;   mu2 = 1.;

  cL.inverse = true;
  cL.A = AcNum/Pe;
  cL.D = 1./Pe;

  OhDerived = active_drop_ohnesorge(Re, Ca);
  mobilityScaleRatio = active_drop_velocity_scale_ratio
    (AcNum, GammaSlope, viscosityRatio, true);
  PeMobility = mobilityScaleRatio*Pe;

  if (pid() == 0)
    fprintf (stderr, "CaseNo=%d Re=%g Ca=%g Pe=%g GammaSlope=%g AcNum=%g "
             "viscosityRatio=%g densityRatio=%g Oh=%g mobilityScaleRatio=%g "
             "PeMobility=%g MAXlevel=%d MINlevel=%d L0=%g tmax=%g tsnap=%g "
             "threshold=%g\n", CaseNo, Re, Ca, Pe, GammaSlope, AcNum,
             viscosityRatio, densityRatio, OhDerived, mobilityScaleRatio, PeMobility,
             MAXlevel, MINlevel, L0, tmax, tsnap, movement_threshold);

  if (pid() == 0)
    fprintf(stderr, "geometry=%s wallHalfWidth=%g dropOffset=%g initialDipole=%g\n",
            "pipe", wallHalfWidth,
            dropOffset, initialDipole);
  run();
  return exit_status;
}

/**
## Initialisation

The unit drop is a sphere centred on the axis.
The optional nonnegative concentration dipole breaks fore-aft symmetry.
Solid geometry and cylindrical metrics precede the CLSVOF initialisation.
*/
event init (i = 0) {
  confined_geometry (wallHalfWidth);
  foreach() {
    double radius = sqrt(sq(x) + sq(y - dropOffset));
    d[] = cs[] > 0. ? 1. - radius : min(1. - radius, -Delta);
    u.x[] = 0.0;
    u.y[] = 0.0;
    cL[] = cs[] > 0. && radius > 1. ?
      initialDipole*(1. + sin(2.*pi*x/L0))*exp(1. - radius) : 0.;
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
  adapt_wavelet ({cs, f, u.x, u.y, cL, KAPPA},
                 (double[]){1e-3, FErr, VelErr, VelErr, cErr, KErr},
                 MAXlevel, MINlevel);
  confined_geometry (wallHalfWidth);
  foreach()
    if (cs[] <= 0.) {
      f[] = cL[] = 0.;
      d[] = min(d[], -Delta);
    }
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
  int invalid_field = 0;
  foreach (reduction(+:ke) reduction(+:drop_volume)
           reduction(+:x_moment) reduction(+:y_moment)
           reduction(max:invalid_field)) {
    double ff = clamp (f[], 0., 1.);
    if (cs[] > 0. && (!isfinite(cL[]) || !isfinite(sigmaf[]) || sigmaf[] <= 0.))
      invalid_field = 1;
    ke += 0.5*rho(ff)*(sq(u.x[]) + sq(u.y[]))*dv();
    drop_volume += ff*dv();
    x_moment += ff*periodic_coordinate(x, xcm_last)*dv();
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
  // The radial moment is not a transverse centre-of-mass displacement.
  dist_last = fabs(xcm_last - x0_cm);

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

  if (invalid_field || !isfinite (ke) || (i > 10 && ke >= keLimit)) {
    if (pid() == 0)
      fprintf (stderr, "Invalid species/tension or kinetic energy %.8e at i=%d, t=%g: stopping.\n",
               ke, i, t);
    print_status ("FAILED", i, t);
    exit_status = 3;
    return 1;
  }

  // Keep the interface stencil separated from the solid wall.
  int unresolved_gap = 0;
  double gap_limit = 3.*L0/(1 << MAXlevel);
  foreach (reduction(max:unresolved_gap))
    if (cs[] > 0. && f[] > 1e-6 &&
        wallHalfWidth - fabs(y) - Delta/2. < gap_limit)
      unresolved_gap = 1;
  if (unresolved_gap) {
    if (pid() == 0)
      fprintf(stderr, "Unresolved drop-wall gap.\n");
    print_status("FAILED", i, t);
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
