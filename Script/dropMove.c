/**
 * Simulation of a two-phase droplet system with surfactant effects
 * 
 * This code simulates a droplet with surfactant-modified surface tension
 * using the CLSVOF (Coupled Level Set and Volume of Fluid) method.
 * The system is non-dimensionalized with the following parameters:
 * - Oh (Ohnesorge number): Ratio of viscous to inertial and surface tension forces
 * - Pe (Peclet number): Ratio of advection to diffusion of surfactants
 * - Ca: Lower Ca means more circular drop (surface tension dominates)
 * - AcNum: Constant surfactant flux from the interface of the drop
 */

#define MIN_LEVEL 0
#define MAX_LEVEL 9

#define VelErr 1e-3
#define FErr 1e-3
#define cErr 1e-3
#define KErr 1e-3

#define tsnap 1e-1

#include "navier-stokes/centered.h"
#define FILTERED
#include "two-phase-clsvof.h"
#include "integral.h"
#include "src-local/activity.h"
// #include "curvature.h"

/**
 * Global variables and boundary conditions
 * cL: Surfactant concentration field
 * sigmaf: Surface tension coefficient field
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

scalar * list = NULL;
int ny, nx; 
double Deltay, Deltax;
double dtmax, tmax; 
double Pe;

scalar sigmaf[];

/**
 * Non-dimensional parameters
 */
#define Oh 1e0
#define Ca 0.1
#define AcNum 1e0

/**
 * Main function: Sets up and runs the simulation
 */
int main(int argc, char const *argv[]){
  stokes = true;
  L0 = 10.0;
  origin (-0.5*L0, -0.5*L0);
  Pe = atof(argv[1]);
  N = 1 << MAX_LEVEL;
  init_grid (N);

  d.sigmaf = sigmaf;

  rho1 = 4/sq(Oh); rho2 = 4/sq(Oh);
  tmax = 50.;
  mu1 = 1.0; mu2 = 1.0;

  cL.inverse = true;
  cL.A = AcNum;
  cL.D = 1e0/Pe;

  char comm[160];
  sprintf (comm, "rm -rf intermediate");
  system(comm);
  sprintf (comm, "mkdir -p intermediate");
  system(comm);
  
  run();  

}

/**
 * Initialization event
 * Sets up initial conditions for:
 * - Distance function (d)
 * - Velocity fields (u.x, u.y)
 * - Surfactant concentration (cL)
 * - Surface tension coefficient (sigmaf)
 */
event init (i = 0) {
  // fraction (fphi, sq(rin) - sq(x) - sq(y));
  foreach() {
    d[] = 1. - sqrt (sq(x) + sq(y));
    u.x[] = 0.0;
    u.y[] = 0.0;
    cL[] = c0; //sq(x) + sq(y) > sq(1.) ? (sq(x) + sq(y) > sq(rout) ? c0 : c0 - cL.A*log(sqrt(sq(x) + sq(y))/rout)) : 0.;
    sigmaf[] = 1/Ca + 4*cL[];
  }
}

/**
 * Properties update event
 * Updates the surface tension coefficient based on surfactant concentration
 */
event properties(i++){
  foreach(){
    sigmaf[] = 1/Ca + 4*cL[];
  }
}

scalar KAPPA[];
event adapt(i++){  

  foreach()
    KAPPA[] = distance_curvature(point, d);
  
  adapt_wavelet({f, u.x, u.y, cL, KAPPA}, (double[]){FErr, VelErr, VelErr, cErr, KErr}, MAX_LEVEL, MIN_LEVEL);
}

/**
 * Output event
 * Saves simulation snapshots at regular intervals
 */
event outputs (t = 0.; t += tsnap; t <= tmax) { 
  char dumpFile[160];
  sprintf (dumpFile, "intermediate/snapshot-%5.4f", t);
  dump (file = dumpFile);
}
scalar cTest[];

/**
 * Logging event
 * Computes and logs kinetic energy of the system
 * Also performs assertions to check simulation stability
 */
event logWriting (i++) {
  
  double ke = 0.;
  foreach(reduction(+:ke)){
    ke += 0.5*rho(f[])*(sq(u.x[])+sq(u.y[]))*sq(Delta);
  }

  double sumv1 = 0.;
  double sumv2 = 0.;
  double sumf = 0.;
  double xcm1 = 0.;
  double ycm2 = 0.;
  double dist = 0.;

  foreach(reduction(+:sumv1), reduction(+:sumv2), reduction(+:sumf)) {
    sumv1 += clamp(f[], 0., 1.)*x;
    sumv2 += clamp(f[], 0., 1.)*y;

    sumf += clamp(f[], 0., 1.);    
  }

  xcm1 = sumv1/sumf;
  ycm2 = sumv2/sumf;

  dist = sqrt(sq(xcm1) + sq(ycm2));

  static FILE * fp;
  if (pid() == 0){
    if (i == 0){
      fprintf (ferr, "i t ke\n");
      fp = fopen ("log.dat", "w");
      fprintf (fp, "i t ke\n");
      fclose (fp);
    }
    fprintf (ferr, "%d %g %5.5e %5.5e\n", i, t, ke, dist);
    fp = fopen ("log.dat", "a");
    fprintf (fp, "%d %g %5.5e %5.5e\n", i, t, ke, dist);
    fclose (fp);
  }

  if (i > 10 && dist >= 1.0) {
    if (pid() == 0) {
      fprintf(stdout, "STATUS MOVED\n");
      fflush(stdout);
    }
    return 1; 
  }

  if (i > 10){
    assert(ke < 1e3);
  }
}

event end (t = tmax) {
  if (pid() == 0) {
    fprintf(stdout, "STATUS NOT_MOVED\n");
    fflush(stdout);
  }
  return 1;
}