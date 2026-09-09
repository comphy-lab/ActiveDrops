# ActiveDrops: repository guide

Basilisk C implementation of a planar chemically active drop and of a
bracketed search for its finite-time onset of self-propulsion in the Péclet
number. Read this file before changing code; `README.md` is the entry point
for building and running.

## Layout

- `simulationCases/dropMove.c`: the planar reference entry point. Reads a
  `key=value` parameter file through `src-local/params.h` and writes
  `intermediate/snapshot-<t>`, `log.dat`, one `STATUS` line and one `SUMMARY`
  line. Generated case folders `simulationCases/c<CaseNo>/` and
  `simulationCases/pescan-<tag>/` are ignored by Git.
- `simulationCases/dropMove-embed-pipe.c` and `dropMove-embed-channel.c`:
  confined axisymmetric and planar cases, sharing `src-local/dropMove-embed.h`.
  Use their matching `embed-*.params` files. Inert embedded walls impose
  no slip and zero species flux; no wetting/contact model is provided.
- `src-local/`: `activity.h` (interfacial chemical source and species
  transport), `parse_params.h` and `params.h` (the single runtime-parameter
  pathway), `two-phase-clsvof-VP.h` (experimental viscoplastic variant, not
  included by the driver).
- `postProcess/`: snapshot readers (`get*.c`, compiled with `qcc`) and the
  plotting scripts `contour.py` and `vectors.py`.
- `testCases/`: software tests only (see evidence classes below).
- Root: `runSimulation.sh`, `runParameterSweep.sh`, `PeScan.py`,
  `default.params`, `sweep.params`.
- `basilisk/` and `.project_config`: the project-local pinned Basilisk
  installed by the CoMPhy ref-locked script. Never committed; the runners
  source `.project_config`.

## Conventions

- One parameter pathway. Runtime values enter only through the `.params`
  file passed as `argv[1]`; do not add `key=value` command-line parsing or a
  second parser. New parameters get a default in the relevant driver, a line in
  `default.params` and a row in the driver's header table.
- Cases run in their own directory. `runSimulation.sh` copies the source and
  parameter file into `simulationCases/c<CaseNo>/`, compiles with
  `-I../../src-local` and executes there. `CaseNo >= 1000`.
- Sweeps go through `runParameterSweep.sh` (`SWEEP_*` Cartesian product,
  deterministic `CaseNo`, `--dry-run` first). Onset scans go through
  `PeScan.py`, which calls `runSimulation.sh` per sample.
- Diagnostics use volume weights. Centroid, drop volume and kinetic energy
  are reduced with `f*dv()`; a cell-count weighted sum is a defect on an
  adaptive mesh. Displacement is measured from the initial centroid.
- Classification is a convention, not a critical value. `MOVED` means the
  centroid displacement exceeded `threshold` before `tmax` at the given
  `MAXlevel`. Any quoted transition interval must carry those three values;
  `PeScan.py` records them in `results.json`.
- Basilisk compatibility. `FILTERED` is defined with a value, reductions in
  `foreach` are space separated, and the removed `dirty` attribute is not
  set. Before asserting Basilisk API behaviour, check the source; do not
  guess.
- Literate C. Documentation lives in `/** ... */` Markdown blocks and is
  rendered by `.github/scripts/generate_docs.py` into `.github/docs/`. Keep
  the file header and section headings when editing.

## Evidence classes

- Software tests (`testCases/`): `test_pescan.py` checks the search logic
  against synthetic classifiers (monotone threshold, all moving, all
  stationary, non-monotone window, numerical failure, run cap).
  `centroid-check.c` checks the centroid diagnostic against the exact
  centroid of a circle on a deliberately asymmetric adaptive mesh. Both
  establish implementation contracts only.
- No verification or validation case exists yet. A single-drop run that
  reports `MOVED` demonstrates that the code runs and that the instability
  develops at that resolution; it does not verify convergence to the
  continuum problem or validate against independent data. Do not describe
  scan output as either.

## Workflow

```sh
bash testCases/run-tests.sh            # software tests (needs qcc for the C check)
bash testCases/run-embed-tests.sh      # embedded geometry and species wall checks
bash runSimulation.sh default.params   # one case
bash runParameterSweep.sh sweep.params --dry-run
python3 PeScan.py 1.0 0.5 --tmax 50 --max-level 9 --tag level9
bash .github/scripts/build.sh          # regenerate the docs site
```

Commit the regenerated `.github/docs/` together with the source change that
motivated it. Do not commit `basilisk/`, `.project_config`, `CLAUDE.md`,
case folders or compiled binaries.
