#!/usr/bin/env python3
"""Bracketed search for the finite-time onset of self-propulsion in Pe.

Each sample is one case run through ``runSimulation.sh`` with a generated
parameter file (the base file plus ``CaseNo``, ``Pe``, ``tmax``,
``MAXlevel``, ``tsnap`` and ``threshold``), so every scan case lives in
``simulationCases/c<CaseNo>/`` exactly like a hand-run case. The driver
classifies a run as ``MOVED`` (the drop centroid moved by more than
``threshold`` before ``tmax``), ``NOT_MOVED`` (it did not within the
observation window) or ``FAILED`` (the run stopped on a numerical failure).

The search first establishes a bracket, one ``NOT_MOVED`` and one ``MOVED``
endpoint, then bisects it, then takes verification samples outside the
bracket. It reports the transition as an interval ``[pe_lo, pe_hi]`` with
``pe_lo`` stationary and ``pe_hi`` moving, never as a rounded point value.
If no bracket exists inside ``[pe_min, pe_max]`` the result is explicitly
``undetermined``.

Every run is recorded (Pe, CaseNo, status, final time, final displacement
and the driver's complete scale-aware SUMMARY line) in a JSON results file so
that the classification convention, resolution, material coupling and
observation horizon travel with the number. The scan varies `Pe`; interpret it
as the mobility-based `PeMobility` only when the fixed base-file scale ratio is
one.

The classification is a finite-time, finite-displacement convention: near
onset slow growth can be censored by ``tmax``. Refine ``tmax``,
``threshold`` and ``MAXlevel`` near the transition and inspect each case's
``log.dat`` for the displacement history before quoting a value.
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import math
import re
import subprocess
import sys
from pathlib import Path
from typing import Callable, Dict, List, Optional

MOVED = "MOVED"
NOT_MOVED = "NOT_MOVED"
FAILED = "FAILED"
STATUSES = (MOVED, NOT_MOVED, FAILED)

REPO_ROOT = Path(__file__).resolve().parent


@dataclasses.dataclass
class Sample:
    pe: float
    status: str
    summary: Optional[Dict[str, str]] = None
    case_no: Optional[int] = None


@dataclasses.dataclass
class ScanResult:
    outcome: str                   # "bracketed" | "undetermined" | "failed" | "nonmonotone"
    pe_lo: Optional[float]         # largest Pe classified NOT_MOVED below the bracket
    pe_hi: Optional[float]         # smallest Pe classified MOVED above the bracket
    reason: str
    samples: List[Sample]
    tolerance_reached: bool = False

    @property
    def width(self) -> Optional[float]:
        if self.pe_lo is None or self.pe_hi is None:
            return None
        return self.pe_hi - self.pe_lo

    def to_dict(self) -> dict:
        return {
            "outcome": self.outcome,
            "pe_lo": self.pe_lo,
            "pe_hi": self.pe_hi,
            "width": self.width,
            "tolerance_reached": self.tolerance_reached,
            "reason": self.reason,
            "samples": [dataclasses.asdict(s) for s in self.samples],
        }


Classifier = Callable[[float], Sample]


def parse_driver_output(out: str) -> Sample:
    """Extract the single STATUS line and the SUMMARY key=value line."""
    status = None
    summary = None
    for raw in out.splitlines():
        line = raw.strip()
        if line.startswith("STATUS "):
            value = line[len("STATUS "):].strip()
            if value not in STATUSES:
                raise RuntimeError(f"Unknown STATUS value '{value}'")
            if status is not None:
                raise RuntimeError("Driver printed more than one STATUS line")
            status = value
        elif line.startswith("SUMMARY "):
            summary = {}
            for token in line[len("SUMMARY "):].split():
                if "=" in token:
                    key, value = token.split("=", 1)
                    summary[key] = value
    if status is None:
        tail = "\n".join(out.splitlines()[-20:])
        raise RuntimeError(f"Did not find STATUS line in output.\n--- tail ---\n{tail}")
    pe = math.nan
    if summary and "Pe" in summary:
        try:
            pe = float(summary["Pe"])
        except ValueError as exc:
            raise RuntimeError(f"Unparsable SUMMARY Pe value '{summary['Pe']}'") from exc
    return Sample(pe=pe, status=status, summary=summary)


def write_case_params(base: Path, target: Path, overrides: Dict[str, object]) -> None:
    """Copy ``base`` to ``target`` with ``key=value`` overrides applied.

    Existing keys are replaced in place; missing keys are appended. Comments
    and unrelated lines are preserved.
    """
    lines = base.read_text().splitlines()
    remaining = dict(overrides)
    out: List[str] = []
    for line in lines:
        m = re.match(r"^\s*([A-Za-z0-9_]+)\s*=", line)
        if m and m.group(1) in remaining:
            key = m.group(1)
            out.append(f"{key}={remaining.pop(key)}")
        else:
            out.append(line)
    for key, value in remaining.items():
        out.append(f"{key}={value}")
    target.write_text("\n".join(out) + "\n")


def make_runner_classifier(base_params: Path, case_start: int, tmax: float, tsnap: float,
                           max_level: int, threshold: float, exec_name: str,
                           threads: int, scan_dir: Path, verbose: bool = True) -> Classifier:
    """Classifier that runs one case per Pe through ``runSimulation.sh``.

    Case numbers start at ``case_start`` and increase by one per sample; the
    generated parameter files are kept in ``scan_dir`` next to the results.
    """
    state = {"case_no": case_start}
    runner = REPO_ROOT / "runSimulation.sh"

    def classify(pe: float) -> Sample:
        case_no = state["case_no"]
        state["case_no"] += 1
        params = scan_dir / f"case_{case_no}.params"
        write_case_params(base_params, params, {
            "CaseNo": case_no, "Pe": f"{pe:.10g}", "tmax": f"{tmax:.10g}",
            "tsnap": f"{tsnap:.10g}", "MAXlevel": max_level,
            "threshold": f"{threshold:.10g}",
        })
        cmd = ["bash", str(runner), str(params), "--exec", exec_name,
               "--threads", str(threads)]
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                              universal_newlines=True, cwd=REPO_ROOT)
        (scan_dir / f"case_{case_no}.stdout.txt").write_text(proc.stdout)
        (scan_dir / f"case_{case_no}.stderr.txt").write_text(proc.stderr)
        try:
            sample = parse_driver_output(proc.stdout)
        except RuntimeError as exc:
            # A crashed or unparsable run is not classifiable. Record it as FAILED so
            # the search stops with every earlier sample intact in the results file.
            print(f"Pe={pe:.6g} (case {case_no}) -> unparsable driver output "
                  f"(exit {proc.returncode}: {exc}); recorded as {FAILED}, see "
                  f"{scan_dir / f'case_{case_no}.stderr.txt'}", flush=True)
            return Sample(pe=pe, status=FAILED, case_no=case_no,
                          summary={"error": f"exit={proc.returncode}"})
        sample.pe = pe
        sample.case_no = case_no
        if verbose:
            t_end = sample.summary.get("t_end", "?") if sample.summary else "?"
            d_end = sample.summary.get("dist_end", "?") if sample.summary else "?"
            print(f"Pe={pe:.6g} (case {case_no}) -> {sample.status} "
                  f"(t_end={t_end}, dist_end={d_end})", flush=True)
        return sample
    return classify


def find_transition(classify: Classifier, pe_start: float, step: float, tol: float,
                    pe_min: float, pe_max: float, max_runs: int,
                    growth: float = 2.0, verify: int = 1) -> ScanResult:
    """Bracket first, then bisect, then verify.

    ``classify`` maps Pe to a Sample. Bisection presumes a monotone response
    (stationary below onset, moving above) and cannot by itself detect a
    violation, because every sample inside a bracket is consistent with it.
    After the bracket has converged, ``verify`` additional samples are taken
    on each side, at distances ``step, 2*step, ...`` outside the bracket, and
    checked against every earlier sample. A stationary sample above a moving
    one ends the search with outcome ``nonmonotone`` and the full sample
    table. This is a bounded check, not a proof of monotonicity.
    """
    if not (pe_min < pe_max):
        raise ValueError("pe_min must be smaller than pe_max")
    if not (pe_min <= pe_start <= pe_max):
        raise ValueError("pe_start must lie inside [pe_min, pe_max]")
    if step <= 0 or tol <= 0 or max_runs < 1:
        raise ValueError("step, tol and max_runs must be positive")

    samples: List[Sample] = []
    stationary: Dict[float, Sample] = {}
    moving: Dict[float, Sample] = {}

    def record(pe: float) -> Sample:
        s = classify(pe)
        s.pe = pe
        samples.append(s)
        if s.status == NOT_MOVED:
            stationary[pe] = s
        elif s.status == MOVED:
            moving[pe] = s
        return s

    def consistent() -> bool:
        if not stationary or not moving:
            return True
        return max(stationary) < min(moving)

    def result(outcome: str, reason: str, tol_ok: bool = False) -> ScanResult:
        lo = max(stationary) if stationary else None
        hi = min(moving) if moving else None
        return ScanResult(outcome=outcome, pe_lo=lo, pe_hi=hi, reason=reason,
                          samples=samples, tolerance_reached=tol_ok)

    def finalise_failure(s: Sample) -> ScanResult:
        return result("failed",
                      f"run at Pe={s.pe:.6g} reported STATUS FAILED; the response is "
                      "not classifiable there, so no bracket is claimed")

    # ---- Phase 1: establish a bracket ----------------------------------
    s = record(pe_start)
    if s.status == FAILED:
        return finalise_failure(s)
    pe = pe_start
    current_step = step
    hit_min = hit_max = False
    while not (stationary and moving):
        if len(samples) >= max_runs:
            return result("undetermined",
                          f"no bracket after {len(samples)} runs (max_runs reached); "
                          f"sampled range [{min(x.pe for x in samples):.6g}, "
                          f"{max(x.pe for x in samples):.6g}]")
        direction = -1.0 if moving else 1.0      # moving -> search lower Pe
        candidate = pe + direction*current_step
        if candidate <= pe_min:
            if hit_min:
                return result("undetermined",
                              f"every sampled Pe down to pe_min={pe_min:g} moved; the onset, "
                              "if any, lies below the search floor")
            candidate, hit_min = pe_min, True
        if candidate >= pe_max:
            if hit_max:
                return result("undetermined",
                              f"every sampled Pe up to pe_max={pe_max:g} stayed stationary; "
                              "the onset, if any, lies above the search ceiling or beyond tmax")
            candidate, hit_max = pe_max, True
        s = record(candidate)
        if s.status == FAILED:
            return finalise_failure(s)
        pe = candidate
        current_step *= growth

    # ---- Phase 2: bisect the bracket -----------------------------------
    tol_ok = False
    while True:
        lo, hi = max(stationary), min(moving)
        if hi - lo <= tol:
            tol_ok = True
            break
        if len(samples) >= max_runs:
            return result("bracketed",
                          f"max_runs={max_runs} reached with bracket width {hi - lo:.6g} > tol={tol:g}")
        mid = 0.5*(lo + hi)
        s = record(mid)
        if s.status == FAILED:
            return finalise_failure(s)

    # ---- Phase 3: bounded monotonicity verification outside the bracket --
    sampled = {x.pe for x in samples}
    verified = 0
    for k in range(1, verify + 1):
        for candidate in (max(stationary) - k*step, min(moving) + k*step):
            candidate = min(max(candidate, pe_min), pe_max)
            if candidate in sampled or len(samples) >= max_runs:
                continue
            sampled.add(candidate)
            verified += 1
            s = record(candidate)
            if s.status == FAILED:
                return finalise_failure(s)
            if not consistent():
                return result("nonmonotone",
                              f"verification sample at Pe={candidate:.6g} is inconsistent with the "
                              "bracket: a stationary sample lies above a moving one, so the "
                              "response is not monotone under this classification convention")
    lo, hi = max(stationary), min(moving)
    return result("bracketed",
                  f"transition bracketed to width {hi - lo:.6g} <= tol={tol:g}; "
                  f"{verified} verification samples consistent",
                  tol_ok=tol_ok)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("pe_start", nargs="?", type=float, default=1.0,
                   help="starting Pe (default 1.0)")
    p.add_argument("step", nargs="?", type=float, default=0.5,
                   help="initial bracketing step (default 0.5); doubles each expansion")
    p.add_argument("--tol", type=float, default=0.005,
                   help="stop when the bracket width is <= tol (default 0.005)")
    p.add_argument("--pe-min", type=float, default=0.001)
    p.add_argument("--pe-max", type=float, default=100.0)
    p.add_argument("--max-runs", type=int, default=60,
                   help="hard cap on the number of simulations (default 60)")
    p.add_argument("--verify", type=int, default=1,
                   help="extra samples per side outside the converged bracket used to check "
                        "monotonicity (default 1)")
    p.add_argument("--base", default="default.params",
                   help="base parameter file copied into every case (default default.params)")
    p.add_argument("--case-start", type=int, default=2000,
                   help="CaseNo of the first scan case; increments by one per run (default 2000)")
    p.add_argument("--tmax", type=float, default=50.0,
                   help="observation horizon written into every case (default 50)")
    p.add_argument("--tsnap", type=float, default=1.0,
                   help="snapshot interval written into every case (default 1)")
    p.add_argument("--max-level", type=int, default=9,
                   help="MAXlevel written into every case (default 9)")
    p.add_argument("--threshold", type=float, default=1.0,
                   help="centroid displacement, in drop radii, that counts as MOVED (default 1)")
    p.add_argument("--exec", default="dropMove.c",
                   help="simulation source in simulationCases/ (default dropMove.c)")
    p.add_argument("--threads", type=int, default=1,
                   help="OpenMP threads per case passed to runSimulation.sh (default 1)")
    p.add_argument("--tag", default="pescan",
                   help="name of the scan; results go to simulationCases/pescan-<tag>/ (default pescan)")
    return p


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    if args.case_start < 1000:
        print("--case-start must be >= 1000 (runner convention)", file=sys.stderr)
        return 2
    if args.threads < 1:
        print("--threads must be a positive integer", file=sys.stderr)
        return 2
    base = Path(args.base)
    if not base.is_absolute():
        base = REPO_ROOT / base
    if not base.is_file():
        print(f"base parameter file not found: {base}", file=sys.stderr)
        return 2
    for line in base.read_text().splitlines():
        if re.match(r"^\s*Oh\s*=", line):
            print("parameter 'Oh' is retired; supply Re and Ca in the base file",
                  file=sys.stderr)
            return 2
    scan_dir = REPO_ROOT / "simulationCases" / f"pescan-{args.tag}"
    scan_dir.mkdir(parents=True, exist_ok=True)
    results_path = scan_dir / "results.json"

    classify = make_runner_classifier(base, args.case_start, args.tmax, args.tsnap,
                                      args.max_level, args.threshold, args.exec,
                                      args.threads, scan_dir)
    result = find_transition(classify, args.pe_start, args.step, args.tol,
                             args.pe_min, args.pe_max, args.max_runs,
                             verify=args.verify)

    try:
        base_label = str(base.relative_to(REPO_ROOT))
    except ValueError:
        base_label = str(base)
    payload = result.to_dict()
    payload["convention"] = {
        "tmax": args.tmax, "threshold": args.threshold, "MAXlevel": args.max_level,
        "tsnap": args.tsnap, "pe_min": args.pe_min, "pe_max": args.pe_max,
        "tol": args.tol, "exec": args.exec, "base_params": base_label,
        "case_start": args.case_start,
    }
    results_path.write_text(json.dumps(payload, indent=2))

    print()
    print(f"Outcome: {result.outcome}")
    print(f"Reason:  {result.reason}")
    if result.pe_lo is not None:
        print(f"Largest stationary Pe: {result.pe_lo:.6g}")
    if result.pe_hi is not None:
        print(f"Smallest moving Pe:    {result.pe_hi:.6g}")
    if result.outcome == "bracketed":
        print(f"Finite-time transition interval (tmax={args.tmax:g}, threshold={args.threshold:g}, "
              f"MAXlevel={args.max_level}): [{result.pe_lo:.6g}, {result.pe_hi:.6g}], "
              f"width {result.width:.6g}")
    print(f"Runs: {len(result.samples)}; results written to {results_path}")
    return 0 if result.outcome == "bracketed" and result.tolerance_reached else 1


if __name__ == "__main__":
    sys.exit(main())
