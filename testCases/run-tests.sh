#!/bin/bash
# run-tests.sh
#
# Run the ActiveDrops software tests from the repository root:
#   1. Python contracts for PeScan, periodic centroids and actual case
#      boundary/solid setup (the case probes require pinned qcc);
#   2. the adaptive-mesh centroid check (needs qcc; skipped with a warning
#      when qcc is unavailable).
#
# Both are software tests in the code-master sense: they establish
# implementation contracts, not convergence or agreement with independent data.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

if [[ -f .project_config ]]; then
  # shellcheck disable=SC1091
  source .project_config
fi

STATUS=0

echo "=== Python software contracts"
if python3 -m unittest discover -s testCases -p 'test_*.py'; then
  echo "PASS: Python software contracts"
else
  echo "FAIL: Python software contracts" >&2
  STATUS=1
fi

echo ""
echo "=== Adaptive-mesh centroid check"
if command -v qcc >/dev/null 2>&1; then
  (
    cd testCases
    qcc -O2 -Wall -disable-dimensions centroid-check.c -o centroid-check -lm
    ./centroid-check
  ) || STATUS=1
else
  echo "WARNING: qcc not found; skipping testCases/centroid-check.c" >&2
fi

exit "$STATUS"
