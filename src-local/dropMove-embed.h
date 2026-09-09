/**
# Embedded confinement driver for a chemically active drop

Shared implementation for `dropMove-embed-pipe.c` and
`dropMove-embed-channel.c`. Uses the baseline CLSVOF/integral surface-tension
model, including its convention $\sigma = 1/Ca + 4c_L$ and source
$A|\nabla f|$. Positive `AcNum` produces species; it does not consume fuel.

The embedded sidewalls impose no slip and $\partial_n c_L=0$. Activity
acts only at the drop interface and adds no solid-wall chemistry. Axial
endcaps impose no slip and `cL=0`; end-distance sensitivity is therefore
part of any confinement study. No contact-angle or wetting model is supplied.
The run stops when a liquid-containing cell enters a three-finest-cell
wall/endcap band.

## Parameters

The baseline parameter-file pathway and output convention are retained.
Additional parameters (lengths in initial drop radii):

| Key | Default | Meaning |
|---|---|---|
| `wallHalfWidth` | 2.53 | Pipe radius or planar channel half-width |
| `dropOffset` | 0 pipe; 0.5 channel | Initial transverse displacement; must be zero in the pipe |
| `initialDipole` | 1e-3 | Nonnegative axial concentration seed; zero preserves fore-aft symmetry |
| `threshold` | 0 | Displacement stop disabled unless supplied by the parameter file |

In the pipe, volume and energy diagnostics use the cylindrical metric but
omit the common factor $2\pi$; `ycm` is the mean radius, not a transverse
centre of mass. The motion classifier uses axial displacement only.

## Author
Vatsal Sanjay (vatsal.sanjay@comphy-lab.org)
Computational Multiphase Physics (CoMPhy) Lab, Durham University
*/

#include <sys/stat.h>
#include <errno.h>

#include "embed.h"
#if ACTIVE_DROP_PIPE
# include "axi.h"
#endif
#include "navier-stokes/centered.h"
#define FILTERED 1
#include "two-phase-clsvof.h"
#include "integral.h"
#include "activity.h"
#include "params.h"
#include "embed-channel-geometry.h"

/**
## Fields and boundary conditions

`cL` is the species concentration in the outer phase (`cL.inverse = true`),
`sigmaf` the surface tension coefficient and `KAPPA` the distance-function
curvature used for adaptation. The embedded sidewalls are inert and no-slip. The axial endcaps are
no-slip concentration reservoirs with `cL = 0`. The pipe axis is a symmetry
boundary, not a solid wall.
*/
scalar cL[], * stracers = {cL};
scalar sigmaf[];
scalar KAPPA[];

cL[embed] = neumann(0.);
cL[top] = neumann(0.);
cL[right] = dirichlet(0.);
cL[left] = dirichlet(0.);
cL[bottom] = neumann(0.);

u.t[top] = dirichlet(0.);
u.t[right] = dirichlet(0.);
u.t[left] = dirichlet(0.);
#if !ACTIVE_DROP_PIPE
u.t[bottom] = dirichlet(0.);
#endif
u.n[embed] = dirichlet(0.);
u.t[embed] = dirichlet(0.);
d[embed] = neumann(0.);
f[embed] = neumann(0.);

/**
## Runtime parameters
*/
int CaseNo = 1000;
double Pe = 1.6;
int MAXlevel = 8, MINlevel = 0;
double Oh = 1., Ca = 0.1, AcNum = 1.;
double tmax = 50., tsnap = 0.1;
double movement_threshold = 0.;
double wallHalfWidth = 2.53, initialDipole = 1e-3;
#if ACTIVE_DROP_PIPE
double dropOffset = 0.;
#else
double dropOffset = 0.5;
#endif
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
  wallHalfWidth = param_double ("wallHalfWidth", wallHalfWidth);
  dropOffset = param_double ("dropOffset", dropOffset);
  initialDipole = param_double ("initialDipole", initialDipole);
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

  double dx = L0_param/(1 << MAXlevel);
  // A face-aligned wall has no cut-cell fragment for embedded no-slip flux.
#if ACTIVE_DROP_PIPE
  double wall_index = wallHalfWidth/dx;
#else
  double wall_index = (wallHalfWidth + L0_param/2.)/dx;
#endif
  if (isfinite(wall_index) && fabs(wall_index - round(wall_index)) < 1e-8) {
    fprintf(stderr, "wallHalfWidth places a wall on a finest-grid face. "
            "Choose a nonaligned value (default 2.53 for L0=10).\n");
    return 2;
  }
  if (!isfinite(wallHalfWidth) || !isfinite(dropOffset) ||
      !isfinite(initialDipole) || initialDipole < 0. ||
      !isfinite(AcNum) || AcNum < 0. || !isfinite(L0_param) ||
      !(wallHalfWidth > 1. + fabs(dropOffset) + 4.*dx) ||
      !(L0_param/2. > 1. + 4.*dx) ||
#if ACTIVE_DROP_PIPE
      dropOffset != 0. || !(wallHalfWidth < L0_param - 2.*dx) ||
#else
      !(wallHalfWidth < L0_param/2. - 2.*dx) ||
#endif
      !(FErr > 0.) || !(VelErr > 0.) || !(cErr > 0.) || !(KErr > 0.) ||
      !isfinite(Pe) || !isfinite(Oh) || !isfinite(Ca) ||
      !isfinite(tmax) || !isfinite(tsnap) || !isfinite(keLimit) ||
      !isfinite(movement_threshold) || !isfinite(FErr) ||
      !isfinite(VelErr) || !isfinite(cErr) || !isfinite(KErr)) {
    fprintf(stderr, "Invalid embedded geometry, dipole or adaptation tolerance. "
            "Keep the unit drop at least four finest cells from walls/endcaps; "
            "pipe dropOffset must be zero.\n");
    return 2;
  }

  if (pid() == 0 && !ensure_directory ("intermediate")) {
    fprintf (stderr, "Cannot create output directory 'intermediate'.\n");
    return 2;
  }

  stokes = true;
  L0 = L0_param;
#if ACTIVE_DROP_PIPE
  origin (-0.5*L0, 0.);
#else
  origin (-0.5*L0, -0.5*L0);
#endif
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

  if (pid() == 0)
    fprintf(stderr, "geometry=%s wallHalfWidth=%g dropOffset=%g initialDipole=%g\n",
            ACTIVE_DROP_PIPE ? "pipe" : "channel", wallHalfWidth,
            dropOffset, initialDipole);
  run();
  return exit_status;
}

/**
## Initialisation

The unit drop is a sphere in the pipe and a circle in the planar channel.
The optional nonnegative concentration dipole breaks fore-aft symmetry.
Solid geometry and cylindrical metrics must precede the CLSVOF initialisation.
*/
static void wall_geometry (void) {
  confined_geometry (wallHalfWidth);
}

event init (i = 0) {
  wall_geometry();
  foreach() {
    double radius = sqrt(sq(x) + sq(y - dropOffset));
    d[] = 1. - radius;
    u.x[] = 0.0;
    u.y[] = 0.0;
    cL[] = cs[] > 0. && radius > 1. ?
      initialDipole*(1. + x/radius)*exp(1. - radius) : 0.;
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
  adapt_wavelet ({cs, f, u.x, u.y, cL, KAPPA},
                 (double[]){1e-3, FErr, VelErr, VelErr, cErr, KErr},
                 MAXlevel, MINlevel);
  wall_geometry();
  foreach()
    if (cs[] <= 0.)
      f[] = cL[] = 0.;
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
#if ACTIVE_DROP_PIPE
  // The radial moment is not a transverse centre-of-mass displacement.
  dist_last = fabs(xcm_last - x0_cm);
#else
  dist_last = sqrt (sq(xcm_last - x0_cm) + sq(ycm_last - y0_cm));
#endif

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

  // The interfacial stencil must remain separated from solids/endcaps.
  int unresolved_gap = 0;
  double gap_limit = 3.*L0/(1 << MAXlevel);
  foreach (reduction(max:unresolved_gap))
    if (cs[] > 0. && f[] > 1e-6 &&
        (wallHalfWidth - fabs(y) - Delta/2. < gap_limit ||
         L0/2. - fabs(x) - Delta/2. < gap_limit))
      unresolved_gap = 1;
  if (unresolved_gap) {
    if (pid() == 0)
      fprintf(stderr, "Unresolved drop-wall/endcap gap; contact is not modelled.\n");
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
