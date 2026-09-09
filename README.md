# ActiveDrops

Spontaneous symmetry breaking of self-propelled drops.

A planar drop emits a chemical species at its interface. The species
diffuses and is advected in the outer phase. The implemented constitutive
law is `sigma = 1/Ca + 4*cL`: positive concentration increases surface
tension. Above a critical Péclet number the isotropic
state is unstable: a small asymmetry in the concentration field drives a
Marangoni flow that reinforces the asymmetry, and the drop self-propels.
The code integrates this problem in the Stokes limit with
[Basilisk](http://basilisk.fr) using the coupled level-set and
volume-of-fluid (CLSVOF) interface method and the integral formulation of
surface tension, and provides a bracketed search for the finite-time onset
in Pe.

## Layout

```
├── basilisk/ - Project-local pinned Basilisk (ignored; installed by the script below)
├── simulationCases/ - Simulation entry point and generated case folders
│   ├── dropMove.c - Planar reference, periodic in both directions
│   ├── dropMove-embed-pipe.c - Axisymmetric drop in a straight embedded pipe
│   └── dropMove-embed-channel.c - Planar drop between embedded walls
├── src-local/ - Project-specific Basilisk headers and the runtime parameter API
│   ├── activity.h - Interfacial chemical source and species transport
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
│   ├── centroid-check.c - Centroid diagnostic on an asymmetric adaptive mesh
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

### Parameter file keys

`dropMove.c` reads these keys (defaults in brackets):

- `CaseNo` [1000], must be at least 1000 so case folders sort.
- `Pe` [1.6]: Péclet number, species diffusivity `1/Pe`.
- `MAXlevel` [8], `MINlevel` [0]: quadtree refinement bounds.
- `Oh` [1], `Ca` [0.1], `AcNum` [1]: density `4/Oh^2`, clean surface tension `1/Ca`, interfacial chemical flux.
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
impermeable to species. Activity adds the source
`AcNum*|grad(f)|` at the drop interface; it does not create wall flux.
Here `cL` is produced species, not directly a consumed-fuel concentration.
Absorbing, reactive or fixed-fuel walls would require a different chemical
boundary model.
Positive activity therefore accumulates species in these periodic,
impermeable domains; there is no imposed chemical sink.

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
