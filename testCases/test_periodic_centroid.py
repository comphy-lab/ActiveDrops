"""
# Periodic centroid software contracts

Compile the actual coordinate helper and moment/displacement expressions
from each case into a small C harness. Synthetic compact weighted samples
cross periodic seams without running the flow solver. These checks assume
that each sample stays within half a period of the preceding centroid.
"""

import math
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import unittest


CASES = ("dropMove", "dropMove-embed-channel", "dropMove-embed-pipe")
REPO = Path(__file__).resolve().parents[1]


def extract_one(pattern, source):
    """Require one production expression, so layout changes fail visibly."""
    matches = re.findall(pattern, source, re.MULTILINE)
    if len(matches) != 1:
        raise AssertionError(f"Expected one match for {pattern!r}, got {len(matches)}")
    return matches[0]


def harness_source(source):
    """Embed production logic; replace only the mesh-volume accessor."""
    helper = extract_one(
        r"(static double periodic_coordinate\s*\([^)]*\)\s*\{[^}]*\})", source
    )
    statements = [extract_one(rf"^\s*({name}\s*\+=\s*[^;]+;)", source)
                  for name in ("drop_volume", "x_moment", "y_moment")]
    displacement = extract_one(r"^\s*(dist_last\s*=\s*[^;]+;)", source)
    return """
#include <math.h>
#include <stdio.h>
#define sq(a) ((a)*(a))
#define dv() volume
double L0;
""" + helper + """
int main(void) {
  double xcm_last, ycm_last, x0_cm, y0_cm;
  if (scanf("%lf %lf %lf %lf %lf", &L0, &xcm_last, &ycm_last,
            &x0_cm, &y0_cm) != 5) return 2;
  int count;
  while (scanf("%d", &count) == 1) {
    double drop_volume = 0., x_moment = 0., y_moment = 0., dist_last;
    for (int j = 0; j < count; j++) {
      double x, y, ff, volume;
      if (scanf("%lf %lf %lf %lf", &x, &y, &ff, &volume) != 4) return 3;
""" + "\n".join(statements) + """
    }
    xcm_last = x_moment/drop_volume;
    ycm_last = y_moment/drop_volume;
""" + displacement + """
    printf("%.17g %.17g %.17g\\n", xcm_last, ycm_last, dist_last);
  }
  return 0;
}
"""


class PeriodicCentroidTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        compiler = shutil.which("cc")
        if compiler is None:
            raise RuntimeError("A C compiler is required for centroid contracts")
        cls.temporary = tempfile.TemporaryDirectory(prefix="active-drops-centroid-")
        cls.addClassCleanup(cls.temporary.cleanup)
        cls.executables = {}
        for case in CASES:
            source = (REPO / "simulationCases" / f"{case}.c").read_text()
            c_file = Path(cls.temporary.name) / f"{case}.c"
            executable = c_file.with_suffix("")
            c_file.write_text(harness_source(source))
            subprocess.run([compiler, "-std=c99", "-Wall", "-Werror",
                            str(c_file), "-o", str(executable), "-lm"],
                           check=True, capture_output=True, text=True, timeout=30)
            cls.executables[case] = executable

    def check_trajectory(self, case, centres, reference=None):
        period = 10.
        initial = centres[0]
        reference = initial if reference is None else reference
        lines = [f"{period} {reference[0]} {reference[1]} {initial[0]} {initial[1]}"]
        fractions = [1., 0.5, 0.75, 1.]
        volumes = [0.125, 0.25, 0.5, 1.]
        weights = [a*b for a, b in zip(fractions, volumes)]
        offsets_x = [-0.3, -0.1, 0.14, 0.36]
        offsets_y = [0.28, -0.32, 0.1, -0.16]
        mean_x = sum(w*x for w, x in zip(weights, offsets_x))/sum(weights)
        mean_y = sum(w*y for w, y in zip(weights, offsets_y))/sum(weights)
        for centre_x, centre_y in centres:
            lines.append(str(len(weights)))
            for dx, dy, fraction, volume in zip(
                    offsets_x, offsets_y, fractions, volumes):
                x = (centre_x + dx - mean_x + period/2.) % period - period/2.
                y = centre_y + dy - mean_y
                if case == "dropMove":
                    y = (y + period/2.) % period - period/2.
                lines.append(f"{x:.17g} {y:.17g} {fraction:.17g} {volume:.17g}")
        result = subprocess.run([str(self.executables[case])],
                                input="\n".join(lines) + "\n", text=True,
                                capture_output=True, check=True, timeout=10)
        actual = [tuple(map(float, row.split())) for row in result.stdout.splitlines()]
        self.assertEqual(len(actual), len(centres))
        previous = None
        for index, ((x, y, distance), (expected_x, expected_y)) in enumerate(
                zip(actual, centres)):
            with self.subTest(case=case, frame=index):
                self.assertAlmostEqual(x, expected_x, delta=1e-11)
                self.assertAlmostEqual(y, expected_y, delta=1e-11)
                dx, dy = expected_x - initial[0], expected_y - initial[1]
                expected_distance = abs(dx) if case.endswith("pipe") else math.hypot(dx, dy)
                self.assertAlmostEqual(distance, expected_distance, delta=1e-11)
                if previous is not None:
                    step_x = expected_x - centres[index - 1][0]
                    step_y = expected_y - centres[index - 1][1]
                    self.assertAlmostEqual(x - previous[0], step_x, delta=1e-11)
                    self.assertAlmostEqual(y - previous[1], step_y, delta=1e-11)
                previous = (x, y)

    def test_baseline_crosses_both_seams_repeatedly(self):
        for direction in (-1., 1.):
            self.check_trajectory("dropMove", [
                (direction*(4.85 + 0.4*i), -direction*(4.9 + 0.3*i))
                for i in range(151)
            ])

    def test_confined_axial_wraps_preserve_transverse_coordinate(self):
        for case in CASES[1:]:
            for direction in (-1., 1.):
                centres = [(direction*(4.85 + 0.4*i), 1.2 + 0.002*i)
                           for i in range(151)]
                # A distant transverse reference must not wrap a solid-wall axis.
                self.check_trajectory(case, centres, reference=(centres[0][0], 101.2))

    def test_wraps_reverse_direction_without_displacement_jump(self):
        for case in CASES:
            path = list(range(101)) + list(range(99, -102, -1))
            self.check_trajectory(case, [
                (0.4*i, -0.3*i if case == "dropMove" else 1.5) for i in path
            ])


if __name__ == "__main__":
    unittest.main()
