# ActiveDrops

Spontaneous symmetry breaking of self-propelled drops.

A planar drop emits a coarse-grained product at its interface. The product
diffuses and is advected in the outer phase. The constitutive law
`sigma = 1/Ca + GammaSlope*cL` makes positive product raise surface
tension. Above a critical Péclet number the isotropic
state is unstable: a small asymmetry in the concentration field drives a
Marangoni flow that reinforces the asymmetry, and the drop self-propels.
The code integrates the dimensionless two-fluid Navier--Stokes equations with
[Basilisk](http://basilisk.fr) using the coupled level-set and
volume-of-fluid (CLSVOF) interface method and the integral formulation of
surface tension. The default `Re=0.01` is small, but creeping-flow behaviour
requires a Reynolds-number convergence study. The code also provides a
bracketed search for the finite-time onset in `Pe`.

## Layout

```
├── basilisk/ - Project-local pinned Basilisk (ignored; installed by the script below)
├── simulationCases/ - Simulation entry point and generated case folders
│   ├── dropMove.c - Planar reference, periodic in both directions
│   ├── dropMove-embed-pipe.c - Axisymmetric drop in a straight embedded pipe
│   └── dropMove-embed-channel.c - Planar drop between embedded walls
├── src-local/ - Project-specific Basilisk headers and the runtime parameter API
│   ├── activity.h - Interfacial chemical source and species transport
│   ├── active-drop-model.h - Scales and mobility diagnostics
│   ├── parse_params.h - Low-level key/value parser for parameter files
│   ├── params.h - Typed parameter accessors (param_int, param_double, ...)
│   └── two-phase-clsvof-VP.h - Experimental viscoplastic CLSVOF variant (not used by dropMove.c)
├── postProcess/ - Snapshot readers and plotting scripts
│   ├── getCM.c - Volume-weighted drop centroid of a snapshot
│   ├── getData.c - Concentration, deformation-rate norm and speed on a grid
│   ├── getDataSlice.c - Volume fraction and velocity on a grid
│   ├── getFacets.c - Interface segments
│   ├── getVelocity_v2.c - Drop velocity of a snapshot
│   ├── contour.py - Parallel per-snapshot panels (concentration, deformation rate, speed)
│   └── vectors.py - Legacy serial velocity-vector frames
├── testCases/ - Software tests
│   ├── test_pescan.py - Synthetic-classifier tests for PeScan.py
│   ├── test_periodic_centroid.py - Production moments across periodic seams
│   ├── test_case_boundaries.py - Actual case setup and solid cleanup probes
│   ├── test_dimensionless_model.py - Production scale mapping and legacy-input rejection
│   ├── centroid-check.c - Centroid diagnostic on an asymmetric adaptive mesh
│   ├── activity-source-budget.c - Discrete interfacial source budget
│   ├── run-tests.sh - Runs Python contracts and the centroid check
│   └── run-embed-tests.sh - Embedded geometry, species and tracer contracts
├── .github/ - Documentation generator, workflows, issue templates and the generated site
├── runSimulation.sh - Single-case compile/run driver
├── runParameterSweep.sh - Parameter sweep driver
├── PeScan.py - Bracketed onset search in Pe (calls runSimulation.sh per sample)
├── default.params - Base runtime parameters
├── sweep.params - Example sweep definition
├── AGENTS.md - Repository conventions and evidence classes
├── LICENSE - GPL-3.0
└── README.md - This file
```

## Requirements

- A C compiler, `make`, `gawk` and `curl` for the Basilisk install.
- Python 3 for `PeScan.py` and the tests (standard library only); the case
  boundary tests also require the pinned Basilisk compiler.
- Python 3 with `numpy`, `pandas` and `matplotlib` for `postProcess/*.py`.

Basilisk is pinned per project. From the repository root:

```sh
curl -fsSL https://raw.githubusercontent.com/comphy-lab/basilisk-C/v2026-08-30/reset_install_basilisk-ref-locked.sh | bash -s -- --ref=v2026-08-30
source .project_config
```

This installs the selected comphy-lab/basilisk-C release into `basilisk/`,
records the ref in `basilisk/.comphy-lock` and writes `.project_config`,
which the runners source. Both paths are ignored by Git. The code was last
built and run against release `v2026-08-30`; it requires a release of
August 2026 or later (`FILTERED` defined with a value, comma-free
reductions, no `dirty` attribute).

## Single case

```sh
bash runSimulation.sh default.params
```

The runner creates `simulationCases/c<CaseNo>/`, copies the parameter and
source files there, compiles with `-I../../src-local` and runs
`./dropMove case.params`. Snapshots go to `intermediate/snapshot-<t>` and
diagnostics to `log.dat` with columns `i t ke dist xcm ycm volume`. The run
prints exactly one `STATUS` line (`MOVED`, `NOT_MOVED` or `FAILED`) and one
`SUMMARY key=value` line. Use `--threads N` for an OpenMP build and
`--exec` to select another source in `simulationCases/`.

## Model and dimensional scales

The initial radius $R_0$, chosen velocity $U_0$ and chosen product
concentration $C_*$ are independent scales. Viscosity is scaled by $\mu_o$,
pressure and stress by $\mu_oU_0/R_0$, and surface tension by $\mu_oU_0$.
The inputs map to the solver as follows.

| Input | Definition | Solver coefficient |
|---|---|---|
| `Re` | $\rho_oU_0R_0/\mu_o$ | $\rho_o=Re$ |
| `Ca` | $\mu_oU_0/\gamma_0$ | reference tension $1/Ca$ |
| `Pe` | $U_0R_0/D$ | product diffusivity $1/Pe$ |
| `GammaSlope` | $\gamma_C C_*/(\mu_oU_0)$ | slope in $1/Ca+\mathit{GammaSlope}\,c_L$ |
| `AcNum` | $A_0R_0/(DC_*)$ | PLIC surface-source coefficient `AcNum/Pe` |
| `viscosityRatio` | $\mu_i/\mu_o$ | $\mu_i=\mathit{viscosityRatio}$ |
| `densityRatio` | $\rho_i/\rho_o$ | $\rho_i=Re\,\mathit{densityRatio}$ |

This is the emitted, surface-tension-increasing interpretation of a
solubilizing drop discussed by
[Michelin (2023)](https://doi.org/10.1146/annurev-fluid-120720-012204).
It does not resolve micelle kinetics, finite fuel or solubilization-driven
drop mass loss; CLSVOF approximately conserves the drop phase.

The reference-tension Ohnesorge number is derived as $Oh=\sqrt{Ca/Re}$. `Oh` is
no longer an input, and a parameter file containing it fails explicitly.

`GammaSlope` is a material coupling under these scales and has the same
default, 4, in all geometries. The mobility comparison is reported separately:

$$
\chi=\frac{U_M}{U_0}=\frac{AcNum\,GammaSlope}{G},\qquad Pe_M=\chi Pe,
$$

where $G=2(1+\lambda)$ for a planar circle and $G=2+3\lambda$ for a sphere,
with $\lambda=\mu_i/\mu_o$. At `AcNum=1`, `GammaSlope=4` and
`viscosityRatio=1`, $\chi=1$ for the planar cases and $0.8$ for the pipe.
Changing only the pipe slope to 5 would normalize its spherical mobility
while changing the material coupling under the same reference scales.
The derivation and evidence limits are in
`src-local/active-drop-model.h`; the planar comparison follows
[Li & Koch (2022)](https://doi.org/10.1017/jfm.2022.891).

### Parameter file keys

`dropMove.c` reads these keys (defaults in brackets):

- `CaseNo` [1000], must be at least 1000 so case folders sort.
- `Re` [0.01]: outer-fluid Reynolds number.
- `Pe` [1.6]: emitted-product Péclet number, diffusivity `1/Pe`.
- `MAXlevel` [8], `MINlevel` [0]: quadtree refinement bounds.
- `Ca` [0.1]: reference surface tension `1/Ca`.
- `GammaSlope` [4]: positive surface-tension sensitivity.
- `AcNum` [1]: outward normal-gradient magnitude; source `AcNum/Pe`.
- `viscosityRatio`, `densityRatio` [1, 1]: inner-to-outer property ratios.
- `L0` [10]: square domain size in drop radii; both directions are periodic.
- `tmax` [50], `tsnap` [0.1]: observation horizon and snapshot interval.
- `threshold` [1]: centroid displacement in drop radii classified as `MOVED`; a value of zero or less disables the early stop.
- `FErr`, `VelErr`, `cErr`, `KErr` [1e-3]: wavelet adaptation tolerances.
- `keLimit` [1e3]: kinetic energy classified as `FAILED`.

## Embedded confinement cases

```sh
bash runSimulation.sh embed-pipe.params --exec dropMove-embed-pipe.c
bash runSimulation.sh embed-channel.params --exec dropMove-embed-channel.c
```

`wallHalfWidth` is the pipe radius or channel half-width, in initial drop
radii (default `2.53`). A wall exactly aligned with finest-grid faces is
rejected because it has no cut-cell fragment for the embedded no-slip flux.
The pipe starts a unit sphere on the axis (`y=0`) and cannot represent
transverse migration. The channel starts a planar unit circle at
`(0,dropOffset)`; its default offset is `0.5`. `initialDipole` seeds an axial
concentration asymmetry and may be set to zero. Both example files disable
the displacement early stop with `threshold=0`.

Each confined case is a standalone simulation driver with the same sections
as `dropMove.c`: fields and boundaries, parameters, `main()`, initialisation,
properties, adaptation, output and diagnostics.

The channel has upper/lower no-slip, non-wetting embedded walls. The pipe
has an upper no-slip, non-wetting embedded wall and a bottom symmetry axis.
Both cases are periodic left/right for every field; there are no axial
endcaps or concentration sinks. The baseline `dropMove.c` is periodic in
both directions. Periodic extent controls interaction with periodic images.

Non-wetting excludes the dispersed phase at the solid: `f[embed]=dirichlet(0)`
and `f=0` in full-solid cells, where the CLSVOF distance remains negative.
Wall chemistry is independent: `cL[embed]=neumann(0)` makes the wall
impermeable to species. Activity adds
$(AcNum/Pe)\,\delta_{\Gamma,h}$ on the reconstructed drop interface; it does
not create wall flux. Pure drop cells and full solids receive no source.
Here `cL` is produced species, not directly a consumed-fuel concentration.
Absorbing, reactive or fixed-fuel walls would require a different chemical
boundary model.
Positive activity therefore accumulates species in these periodic,
impermeable domains; there is no imposed chemical sink.

## Version migration

Replace the retired input block

```ini
Oh=1
Ca=0.1
Pe=1.6
AcNum=1
```

with an explicit scale contract, for example

```ini
Re=0.01
Ca=0.1
Pe=1.6
GammaSlope=4
AcNum=1
viscosityRatio=1
densityRatio=1
```

These examples are not numerically equivalent. The earlier code used
`rho=4/Oh^2`, treated `AcNum` as the diffuse source coefficient and fixed the
tension slope at 4. Historical outputs must retain their original code and
parameter interpretation; they cannot be relabelled with the new groups.

## Model limits

- Product concentration is stored through the two-phase transport
  construction while diffusion is weighted toward the outer phase. A sharp
  exterior-flux limit has not been verified.
- The source uses a geometric PLIC interface measure, including the radial
  metric in axisymmetry. Its planar-circle and spherical-area convergence is
  checked as a software benchmark; a sharp exterior-flux verification remains
  open.
- A spatially uniform product changes the absolute tension and effective
  deformability, although only gradients produce Marangoni stress.
- Periodic impermeable domains retain emitted product. Long-time steady states
  need separate evidence or an explicit chemical relaxation mechanism.
- `stokes=false` retains unsteady and convective inertia consistently. A small
  input `Re` does not replace a Reynolds-number convergence study.
- `MOVED` is a finite-time classifier, not the theoretical $Pe_M=4$ threshold
  for an unbounded, nondeforming sphere.
- Embedded cases stop before contact and contain no moving-contact-line model.

Finite exterior fuel, uptake arrest and delayed chemical sensing belong to
the sister `active-drops-with-memory` programme. Its existing delayed-sampling
PR is independent of this model change.

These cases require the pinned Basilisk `v2026-08-30` and reconstruct the
stationary embedded geometry after mesh adaptation. A liquid-containing
cell entering a three-finest-cell wall band
stops with `FAILED`. In the pipe, `volume` and `ke` omit the common `2*pi`
factor; `ycm` is the mean radius, and `dist` measures axial displacement.
This phase-exclusion treatment does not supply a contact-angle or moving
contact-line law.

Periodic centroid coordinates are unwrapped about the previous centroid,
so a seam crossing does not create a displacement jump; logged coordinates
may leave the principal domain. This assumes a compact single drop within
half a period of the tracked centroid. Snapshot readers need the matching
periodic treatment and, for the pipe, cylindrical volume weights.

Software checks for geometry and species wall conditions are available via
`bash testCases/run-embed-tests.sh`. These are implementation checks;
convergence and physical wall attraction or repulsion require separate study.

## Parameter sweep

```sh
bash runParameterSweep.sh sweep.params --dry-run
bash runParameterSweep.sh sweep.params --parallel 2
```

`SWEEP_*` lists in the sweep file form a Cartesian product; each case gets a
deterministic `CaseNo` from `CASE_START`. Set `CASE_END` to assert the
expected number of cases.

## Onset scan

```sh
python3 PeScan.py 1.0 0.5 --tmax 50 --max-level 9 --tol 0.005 --tag level9
```

`PeScan.py` runs one case per sample through `runSimulation.sh` (case
numbers from `--case-start`, default 2000), first bracketing the transition
with one stationary and one moving endpoint, then bisecting until the
bracket is narrower than `--tol`, then taking one verification sample on
each side. Results, including every sampled run and the classification
convention, are written to `simulationCases/pescan-<tag>/results.json`.
Outcomes are:

- `bracketed`: interval `[pe_lo, pe_hi]` with `pe_lo` stationary and `pe_hi` moving;
- `undetermined`: no bracket inside `[--pe-min, --pe-max]` (all runs moved, or none did within `tmax`);
- `failed`: a run stopped on a numerical failure and is not classifiable;
- `nonmonotone`: a verification sample contradicts a single threshold.

The exit status is 0 only for a bracket converged to tolerance.

## Diagnostics and classification convention

The drop centroid, drop volume and kinetic energy are integrated with the
cell volume, `f*dv()`. On the adaptive mesh the cell count is not a measure
of area, so a count-weighted centroid is biased towards refined regions;
`testCases/centroid-check.c` places a circle on a mesh whose right half is
refined two levels deeper and shows a count-weighted centroid error of
about two coarse cells against below one percent of a coarse cell for the
volume-weighted centroid. Displacement is measured from the centroid at the
first step.

`MOVED` means the centroid moved by more than `threshold` before `tmax`.
This is a finite-time, finite-displacement convention: near onset the
displacement grows slowly and can be censored by the observation window, so
a reported interval depends on `tmax`, `threshold` and `MAXlevel`, all of
which are recorded in every `SUMMARY` line and in `results.json`. Near the
transition, lengthen `tmax`, lower `threshold` and increase `MAXlevel`
until the interval stops moving, and inspect the growth of `dist` in
`log.dat` before quoting a value.

## Tests

```sh
bash testCases/run-tests.sh
```

Both tests are software tests: the first exercises the scan logic with
synthetic classifiers (monotone threshold, all moving, all stationary,
non-monotone window, numerical failure, run cap); the second checks the
centroid diagnostic against the exact centroid of a circle on an
asymmetrically refined quadtree and needs `qcc`. Neither is a verification
of convergence nor a validation against independent data; no such case
exists in this repository yet.

## Post-processing

```sh
cd postProcess
python3 contour.py --caseToProcess ../simulationCases/c1000 --tSnap 0.1 --cpus 4
```

`contour.py` compiles `getFacets` and `getData` once with `qcc`, then
renders one frame per snapshot in parallel (`--cpus`, alias `--CPUs`;
`--max-frames` bounds the number of snapshots). `vectors.py` is a legacy
serial script with hard-coded paths kept for reference. The readers compile
individually with

```sh
qcc -O2 -Wall -disable-dimensions getCM.c -o getCM -lm
./getCM ../simulationCases/c1000/intermediate/snapshot-10.0000
```

## Documentation site

```sh
bash .github/scripts/build.sh
bash .github/scripts/deploy.sh
```

The generator renders the literate C and Python sources under `src-local/`,
`simulationCases/`, `postProcess/` and `testCases/` into `.github/docs/`,
which is committed; the GitHub Pages workflow deploys it from `main`.

## Licence

GPL-3.0; see `LICENSE`.
