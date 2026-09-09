#!/usr/bin/env python3
"""Synthetic-classifier tests for PeScan.find_transition.

Run from the repository root:  python3 -m unittest discover -s testCases -p 'test_*.py' -v
No Basilisk installation is required.
"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import PeScan  # noqa: E402
from PeScan import FAILED, MOVED, NOT_MOVED, Sample, find_transition, parse_driver_output  # noqa: E402


def classifier_from(fn):
    def classify(pe):
        return Sample(pe=pe, status=fn(pe))
    return classify


class FindTransitionTests(unittest.TestCase):
    def test_monotone_threshold_is_bracketed_to_tolerance(self):
        pe_c = 1.37
        calls = []

        def fn(pe):
            calls.append(pe)
            return MOVED if pe > pe_c else NOT_MOVED

        r = find_transition(classifier_from(fn), pe_start=1.0, step=0.5, tol=0.005,
                            pe_min=0.001, pe_max=100.0, max_runs=60)
        self.assertEqual(r.outcome, "bracketed")
        self.assertTrue(r.tolerance_reached)
        self.assertLess(r.pe_lo, pe_c)
        self.assertGreaterEqual(r.pe_hi, pe_c)
        self.assertLessEqual(r.width, 0.005)
        self.assertEqual(len(r.samples), len(calls))
        # every sample is consistent with the bracket
        for s in r.samples:
            if s.status == NOT_MOVED:
                self.assertLessEqual(s.pe, r.pe_lo)
            else:
                self.assertGreaterEqual(s.pe, r.pe_hi)

    def test_monotone_threshold_from_above(self):
        pe_c = 0.42
        r = find_transition(classifier_from(lambda pe: MOVED if pe > pe_c else NOT_MOVED),
                            pe_start=5.0, step=0.5, tol=0.01,
                            pe_min=0.001, pe_max=100.0, max_runs=60)
        self.assertEqual(r.outcome, "bracketed")
        self.assertLess(r.pe_lo, pe_c)
        self.assertGreaterEqual(r.pe_hi, pe_c)
        self.assertLessEqual(r.width, 0.01)

    def test_all_moving_is_undetermined_not_a_value(self):
        r = find_transition(classifier_from(lambda pe: MOVED), pe_start=1.0, step=0.5,
                            tol=0.005, pe_min=0.001, pe_max=100.0, max_runs=60)
        self.assertEqual(r.outcome, "undetermined")
        self.assertIsNone(r.pe_lo)
        self.assertIn("pe_min", r.reason)
        self.assertEqual(min(s.pe for s in r.samples), 0.001)

    def test_all_stationary_is_undetermined_not_a_value(self):
        r = find_transition(classifier_from(lambda pe: NOT_MOVED), pe_start=1.0, step=0.5,
                            tol=0.005, pe_min=0.001, pe_max=100.0, max_runs=60)
        self.assertEqual(r.outcome, "undetermined")
        self.assertIsNone(r.pe_hi)
        self.assertIn("pe_max", r.reason)
        self.assertEqual(max(s.pe for s in r.samples), 100.0)

    def test_nonmonotone_response_is_reported(self):
        # moving only inside a window: stationary below 1, moving in (1, 1.4), stationary above
        def fn(pe):
            return MOVED if 1.0 < pe < 1.4 else NOT_MOVED

        r = find_transition(classifier_from(fn), pe_start=1.2, step=0.5, tol=0.005,
                            pe_min=0.001, pe_max=100.0, max_runs=60)
        self.assertEqual(r.outcome, "nonmonotone")
        self.assertIn("verification", r.reason)

    def test_coarse_scan_can_miss_a_window_and_says_so(self):
        # a window narrower than the step is legitimately reported as undetermined,
        # never as a bracketed value
        def fn(pe):
            return MOVED if 1.0 < pe < 3.0 else NOT_MOVED

        r = find_transition(classifier_from(fn), pe_start=0.5, step=4.0, tol=0.005,
                            pe_min=0.001, pe_max=100.0, max_runs=60)
        self.assertEqual(r.outcome, "undetermined")

    def test_failed_run_stops_the_search(self):
        def fn(pe):
            return FAILED if pe > 2.0 else NOT_MOVED

        r = find_transition(classifier_from(fn), pe_start=1.0, step=0.5, tol=0.005,
                            pe_min=0.001, pe_max=100.0, max_runs=60)
        self.assertEqual(r.outcome, "failed")
        self.assertIsNone(r.pe_hi)

    def test_max_runs_caps_bisection_and_reports_width(self):
        pe_c = 1.37
        r = find_transition(classifier_from(lambda pe: MOVED if pe > pe_c else NOT_MOVED),
                            pe_start=1.0, step=0.5, tol=1e-9,
                            pe_min=0.001, pe_max=100.0, max_runs=6)
        self.assertEqual(r.outcome, "bracketed")
        self.assertFalse(r.tolerance_reached)
        self.assertEqual(len(r.samples), 6)
        self.assertLess(r.pe_lo, pe_c)
        self.assertGreaterEqual(r.pe_hi, pe_c)

    def test_invalid_arguments(self):
        with self.assertRaises(ValueError):
            find_transition(classifier_from(lambda pe: MOVED), pe_start=200.0, step=0.5,
                            tol=0.005, pe_min=0.001, pe_max=100.0, max_runs=60)
        with self.assertRaises(ValueError):
            find_transition(classifier_from(lambda pe: MOVED), pe_start=1.0, step=0.0,
                            tol=0.005, pe_min=0.001, pe_max=100.0, max_runs=60)


class ParseDriverOutputTests(unittest.TestCase):
    def test_status_and_summary_are_parsed(self):
        out = ("i t ke dist\n0 0 0 0\nSTATUS MOVED\n"
               "SUMMARY Pe=2.5 max_level=9 tmax=50 threshold=1 t_end=12.3 i_end=800 "
               "dist_end=1.0012e+00 xcm_end=0.7 ycm_end=0.7 status=MOVED\n")
        s = parse_driver_output(out)
        self.assertEqual(s.status, MOVED)
        self.assertEqual(s.pe, 2.5)
        self.assertEqual(s.summary["t_end"], "12.3")

    def test_missing_status_raises(self):
        with self.assertRaises(RuntimeError):
            parse_driver_output("no status here\n")

    def test_duplicate_status_raises(self):
        with self.assertRaises(RuntimeError):
            parse_driver_output("STATUS MOVED\nSTATUS NOT_MOVED\n")

    def test_unknown_status_raises(self):
        with self.assertRaises(RuntimeError):
            parse_driver_output("STATUS MAYBE\n")


if __name__ == "__main__":
    unittest.main()
