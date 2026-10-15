"""
# Standalone case boundary and solid-state contracts

Compile temporary copies of the actual cases and replace only their time
loop entry with a bounded inspection. The production parameter setup and
metric, defaults, initialisation and adaptation events execute; no timestep
or scientific simulation runs. A pinned project-local qcc is required.
"""

import os
from pathlib import Path
import re
import subprocess
import tempfile
import unittest


REPO = Path(__file__).resolve().parents[1]
CASES = ("dropMove", "dropMove-embed-channel", "dropMove-embed-pipe")

CHECK = r"""
static int inspect_solid (const char * stage)
{
#if EMBED
  int solids = 0, bad = 0;
  foreach(reduction(+:solids) reduction(+:bad)) {
    if (cs[] <= 0.) {
      solids++;
      if (f[] != 0. || cL[] != 0. || ActivityFlux[] != 0. || !(d[] < 0.))
        bad++;
    }
  }
  fprintf(stderr, "solid stage=%s cells=%d invalid=%d\n", stage, solids, bad);
  return solids == 0 || bad != 0;
#else
  return 0;
#endif
}

static int check_case()
{
  event("metric");
  event("defaults");
  event("init");
#if EMBED
  bool periodic_y = false;
#else
  bool periodic_y = true;
#endif
  int bad = !Period.x || Period.y != periodic_y;
  for (scalar field in {cL, f, d, p}) {
    for (int side = 0; side < 4; side++) {
      bool expected = side == left || side == right || periodic_y;
      bool periodic_callback = field.boundary[side] == periodic_bc;
      if (periodic_callback != expected) {
        fprintf(stderr, "boundary scalar=%s side=%d periodic=%d expected=%d\n",
                field.name, side, periodic_callback, expected);
        bad++;
      }
    }
  }
#if EMBED
  int fragments = 0;
  foreach(reduction(+:bad) reduction(+:fragments)) {
    if (cs[] > 0. && cs[] < 1.) {
      fragments++;
      for (scalar field in {cL, f, d, u.x, u.y}) {
        bool dirichlet = false;
        double value = field.boundary[embed](point, point, field, &dirichlet);
        bool expected = field.i == f.i || field.i == u.x.i || field.i == u.y.i;
        if (dirichlet != expected || value != 0.)
          bad++;
      }
    }
  }
  if (fragments == 0) bad++;
  fprintf(stderr, "embedded fragments=%d boundary_errors=%d\n", fragments, bad);
#endif
  /* Exercise the production source before inspecting solid support. */
#if EMBED
  dt = 0.01;
  cL.D = 0.;
  foreach_face() uf.x[] = 0.;
  boundary({f, cL, uf});
  event("vof");
  event("tracer_diffusion");
#endif
  bad += inspect_solid("init");
#if EMBED
  foreach()
    if (cs[] <= 0.) {
      d[] = 0.;
      cL[] = 0.375;
    }
#endif
  event("adapt");
  bad += inspect_solid("adapt");
  fprintf(stderr, "%s: case boundaries periodic_x=%d periodic_y=%d\n",
          bad ? "FAIL" : "PASS", Period.x, Period.y);
  return bad ? 1 : 0;
}
"""


class CaseBoundaryTests(unittest.TestCase):
    def test_actual_case_initialisation_and_adaptation(self):
        compiler = REPO / "basilisk" / "src" / "qcc"
        self.assertTrue(compiler.is_file(), "Pinned project-local qcc is required")
        self.assertTrue((REPO / "basilisk" / ".comphy-lock").is_file(),
                        "Missing project-local Basilisk lock")
        environment = os.environ.copy()
        environment["BASILISK"] = str(compiler.parent)
        environment["PATH"] = str(compiler.parent) + os.pathsep + environment.get("PATH", "")
        with tempfile.TemporaryDirectory(prefix="active-drops-boundaries-") as temporary:
            root = Path(temporary)
            for case in CASES:
                with self.subTest(case=case):
                    work = root / case
                    work.mkdir()
                    source = (REPO / "simulationCases" / f"{case}.c").read_text()
                    source, count = re.subn(r"\brun\s*\(\s*\)\s*;",
                                            "return check_case();", source)
                    self.assertEqual(count, 1, "Expected exactly one time-loop entry")
                    source, count = re.subn(r"\bint main\s*\(",
                                            "static int check_case();\nint main (", source)
                    self.assertEqual(count, 1, "Expected exactly one case main")
                    c_file = work / f"{case}.c"
                    c_file.write_text(source + CHECK)
                    executable = work / "check-case"
                    command = [str(compiler), "-O1", "-disable-dimensions",
                               f"-I{REPO / 'src-local'}", c_file.name,
                               "-o", str(executable), "-lm"]
                    compiled = subprocess.run(command, cwd=work, env=environment,
                                              text=True, capture_output=True, timeout=120)
                    self.assertEqual(compiled.returncode, 0, compiled.stdout + compiled.stderr)
                    parameters = work / "case.params"
                    parameters.write_text(
                        "MAXlevel=6\nMINlevel=3\ntmax=1\nL0=10\n"
                        "wallHalfWidth=2.53\nthreshold=0\n"
                    )
                    result = subprocess.run([str(executable), str(parameters)],
                                            cwd=work, env=environment, text=True,
                                            capture_output=True, timeout=30)
                    self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                    self.assertIn("PASS: case boundaries", result.stderr)
                    print(f"{case}:\n{result.stderr}", end="")


if __name__ == "__main__":
    unittest.main()
