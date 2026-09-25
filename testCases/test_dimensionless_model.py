"""
# Dimensionless parameter mapping contracts

Compile temporary copies of the three production drivers and replace only the
time-loop entry with a direct inspection.  The checks establish the mapping
from public dimensionless inputs to Basilisk coefficients and the reported
mobility diagnostic.  They do not advance the coupled equations.
"""

import math
import os
from pathlib import Path
import re
import subprocess
import tempfile
import unittest


REPO = Path(__file__).resolve().parents[1]
CASES = ("dropMove", "dropMove-embed-channel", "dropMove-embed-pipe")

CHECK = r"""
static int check_model (void)
{
#if AXI
  const int spherical = 1;
#else
  const int spherical = 0;
#endif
  const double expected_denominator = spherical ? 8. : 6.;
  const double expected_mobility = 0.75*6./expected_denominator;
  int bad = stokes || fabs(rho2 - 0.08) > 1e-14 ||
    fabs(rho1 - 0.04) > 1e-14 || fabs(mu2 - 1.) > 1e-14 ||
    fabs(mu1 - 2.) > 1e-14 || fabs(cL.D - 0.2) > 1e-14 ||
    fabs(cL.A - 0.15) > 1e-14 ||
    fabs(OhDerived - sqrt(2.5)) > 1e-14 ||
    fabs(mobilityScaleRatio - expected_mobility) > 1e-14 ||
    fabs(PeMobility - 5.*expected_mobility) > 1e-14;
  fprintf(stderr, "%s: model spherical=%d rho=(%g,%g) mu=(%g,%g) "
          "D=%g source=%g Oh=%g mobility=%g PeMobility=%g\n",
          bad ? "FAIL" : "PASS", spherical, rho1, rho2, mu1, mu2,
          cL.D, cL.A, OhDerived, mobilityScaleRatio, PeMobility);
  return bad ? 1 : 0;
}
"""


class DimensionlessModelTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.compiler = REPO / "basilisk" / "src" / "qcc"
        if not cls.compiler.is_file():
            raise RuntimeError("Pinned project-local qcc is required")
        cls.environment = os.environ.copy()
        cls.environment["BASILISK"] = str(cls.compiler.parent)
        cls.environment["PATH"] = (
            str(cls.compiler.parent) + os.pathsep + cls.environment.get("PATH", "")
        )

    def test_runner_preflight_rejects_empty_retired_parameter(self):
        with tempfile.TemporaryDirectory(prefix="active-drops-retired-") as temporary:
            root = Path(temporary)
            base = root / "base.params"
            sweep = root / "sweep.params"
            sweep.write_text(f"BASE_CONFIG={base}\nSWEEP_Pe=1,2\n")
            for assignment in ("Oh=", "  Oh = 1 # retired"):
                base.write_text(assignment + "\n")
                for script, argument, flags in (
                    ("runSimulation.sh", base, []),
                    ("runParameterSweep.sh", sweep, ["--dry-run"]),
                ):
                    with self.subTest(assignment=assignment, script=script):
                        result = subprocess.run(
                            ["bash", str(REPO / script), str(argument), *flags],
                            cwd=root, env=self.environment, text=True,
                            capture_output=True, timeout=30,
                        )
                        self.assertNotEqual(result.returncode, 0)
                        self.assertIn("Parameter 'Oh' is retired", result.stderr)

    def test_actual_driver_mappings_and_retired_oh(self):
        with tempfile.TemporaryDirectory(prefix="active-drops-model-") as temporary:
            root = Path(temporary)
            for case in CASES:
                with self.subTest(case=case):
                    work = root / case
                    work.mkdir()
                    source = (REPO / "simulationCases" / f"{case}.c").read_text()
                    source, count = re.subn(
                        r"\brun\s*\(\s*\)\s*;", "return check_model();", source
                    )
                    self.assertEqual(count, 1, "Expected one time-loop entry")
                    source, count = re.subn(
                        r"\bint main\s*\(", "static int check_model (void);\nint main (", source
                    )
                    self.assertEqual(count, 1, "Expected one case main")
                    c_file = work / f"{case}.c"
                    c_file.write_text(source + CHECK)
                    executable = work / "check-model"
                    compiled = subprocess.run(
                        [str(self.compiler), "-O1", "-disable-dimensions",
                         f"-I{REPO / 'src-local'}", c_file.name,
                         "-o", str(executable), "-lm"],
                        cwd=work, env=self.environment, text=True,
                        capture_output=True, timeout=120,
                    )
                    self.assertEqual(compiled.returncode, 0,
                                     compiled.stdout + compiled.stderr)
                    params = work / "case.params"
                    params.write_text(
                        "Re=0.08\nCa=0.2\nPe=5\nGammaSlope=6\nAcNum=0.75\n"
                        "viscosityRatio=2\ndensityRatio=0.5\n"
                        "MAXlevel=6\nMINlevel=3\ntmax=1\nL0=10\n"
                        "wallHalfWidth=2.53\nthreshold=0\n"
                    )
                    result = subprocess.run(
                        [str(executable), str(params)], cwd=work,
                        env=self.environment, text=True, capture_output=True, timeout=30,
                    )
                    self.assertEqual(result.returncode, 0,
                                     result.stdout + result.stderr)
                    self.assertIn("PASS: model", result.stderr)

                    legacy = work / "legacy.params"
                    legacy.write_text(params.read_text() + "Oh=1\n")
                    rejected = subprocess.run(
                        [str(executable), str(legacy)], cwd=work,
                        env=self.environment, text=True, capture_output=True, timeout=30,
                    )
                    self.assertEqual(rejected.returncode, 2)
                    self.assertIn("Parameter 'Oh' is retired", rejected.stderr)


if __name__ == "__main__":
    unittest.main()
