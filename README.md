# ActiveDrops

Spontaneous symmetry breaking of self-propelled drops.

A planar drop emits a chemical species at its interface. The species
diffuses and is advected in the outer phase and lowers the interfacial
tension where it accumulates. Above a critical Péclet number the
isotropic state is unstable: a small asymmetry in the concentration field
drives a Marangoni flow that reinforces the asymmetry, and the drop
self-propels. The code integrates this problem in the Stokes limit with
[Basilisk](http://basilisk.fr) using the coupled level-set and
volume-of-fluid (CLSVOF) interface method and the integral formulation of
surface tension.

## Contents

| Path | Purpose |
|---|---|
| `dropMove.c` | Single fixed-Pe demonstration run (Pe = 1.6, level 8). |
| `Script/dropMove.c` | Scan driver: takes Pe on the command line and classifies the run as `MOVED`, `NOT_MOVED` or `FAILED`. |
| `Script/PeScan.py` | Bracketed search for the finite-time onset in Pe; reports an interval, never a rounded point. |
| `Script/tests/` | Synthetic-classifier unit tests for the scan logic and a Basilisk centroid check on an asymmetric adaptive mesh. |
| `src-local/activity.h` | Interfacial chemical source and advection-diffusion of the species (shared by both drivers). |
| `src-local/two-phase-clsvof-VP.h` | Experimental viscoplastic CLSVOF variant; not included by the drivers. |
| `postProcess-contour/`, `postProcess-vectors/` | Snapshot readers (`getData`, `getFacets`, `getCM`, `getVelocity_v2`, `getDataSlice`) and plotting scripts. |
| `runCases.sh` | Compiles and runs the fixed-Pe demonstration. |

## Non-dimensional parameters

The drop radius, the outer viscosity and a chemical flux scale set the
units. The drivers fix `Oh = 1`, `Ca = 0.1` and `AcNum = 1` at compile
time; Pe is a runtime argument.

- `Oh`: Ohnesorge number, entering only through the density `rho = 4/Oh^2`.
- `Pe`: Péclet number; the species diffusivity is `D = 1/Pe`.
- `Ca`: the clean-interface surface tension is `1/Ca`; the local coefficient is `1/Ca + 4 c`, with `c` the species concentration.
- `AcNum`: constant chemical flux emitted at the interface.

The domain is a square of side 10 drop radii with no-slip walls and zero
concentration on all boundaries.

## Requirements

Basilisk with `qcc`, a C compiler and the standard maths library. The
drivers compile against Basilisk releases of August 2026 or later
(`FILTERED` is defined with a value, reductions are listed without commas
and the removed `dirty` attribute is no longer set). A ref-locked Basilisk
can be installed next to this repository with the CoMPhy install script:

```sh
curl -sL https://raw.githubusercontent.com/comphy-lab/basilisk-C/main/reset_install_basilisk-ref-locked.sh | bash -s -- --hard
source .project_config
```

The resulting `basilisk/` directory and `.project_config` are ignored by
Git.

## Build and run

Fixed-Pe demonstration:

```sh
qcc -O2 -Wall -disable-dimensions dropMove.c -o dropMove -lm
./dropMove
```

Scan driver, one Pe value:

```sh
cd Script
qcc -O2 -Wall -disable-dimensions dropMove.c -o dropMove -lm
./dropMove 4 tmax=50 tsnap=1 max_level=9 threshold=1 out=runs/pe-4
```

Optional `key=value` arguments override the observation horizon `tmax`,
the snapshot interval `tsnap`, the maximum refinement level `max_level`,
the displacement `threshold` (in drop radii) that counts as motion, and
the output directory `out`. The run prints exactly one `STATUS` line
(`MOVED`, `NOT_MOVED` or `FAILED`) followed by a `SUMMARY` line with the
parameters, the final time and the final displacement. The displacement
history is written to `<out>/log.dat` with columns
`i t ke dist xcm ycm volume`.

Onset scan:

```sh
cd Script
python3 PeScan.py 1.0 0.5 --tmax 50 --max-level 9 --tol 0.005 --out scan
```

`PeScan.py` first brackets the transition with one stationary and one
moving endpoint, then bisects until the bracket is narrower than `--tol`,
and finally takes one verification sample on each side of the bracket to
check monotonicity. The result is written to `scan/results.json` together
with every sampled run. Outcomes are:

- `bracketed`: interval `[pe_lo, pe_hi]` with `pe_lo` stationary and `pe_hi` moving;
- `undetermined`: no bracket inside `[--pe-min, --pe-max]` (all runs moved, or none did within `tmax`);
- `failed`: a run stopped on a numerical failure and is not classifiable;
- `nonmonotone`: a verification sample contradicts a single threshold.

The exit status is 0 only for a bracket converged to tolerance.

## Diagnostics and classification convention

The drop centroid and kinetic energy are integrated with the cell volume,
`f*dv()`. On the adaptive mesh the cell count is not a measure of area, so
a count-weighted centroid is biased towards refined regions; the test in
`Script/tests/centroid-check.c` places a stationary circle on a mesh whose
right half is refined two levels deeper and shows a count-weighted
centroid error of about two coarse cells against below one percent of a
coarse cell for the volume-weighted centroid. Displacement is measured
from the centroid recorded at the first step.

`MOVED` means the centroid moved by more than `threshold` before `tmax`.
This is a finite-time, finite-displacement convention: near onset the
displacement grows slowly and can be censored by the observation window,
so a reported interval depends on `tmax`, `threshold` and `max_level`.
These are recorded in every `SUMMARY` line and in `results.json`. Near the
transition, lengthen `tmax`, lower `threshold` and increase `max_level`
until the interval stops moving, and inspect the growth of `dist` in
`log.dat` before quoting a value.

## Tests

```sh
cd Script
python3 -m unittest tests.test_pescan -v
qcc -O2 -Wall -disable-dimensions tests/centroid-check.c -o centroid-check -lm && ./centroid-check
```

The first exercises the scan logic with synthetic classifiers
(monotone threshold, all moving, all stationary, non-monotone window,
numerical failure, run cap). The second is the adaptive-mesh centroid
check described above and requires Basilisk.

## Post-processing

The readers in `postProcess-contour/` and `postProcess-vectors/` compile
with `qcc` in the same way as the drivers, for example

```sh
cd postProcess-vectors
qcc -O2 -Wall -disable-dimensions getCM.c -o getCM -lm
./getCM ../Script/runs/pe-4/snapshot-10.0000
```

`getCM` prints the volume-weighted centroid and the snapshot time.

## Licence

GPL-3.0; see `LICENSE`.
