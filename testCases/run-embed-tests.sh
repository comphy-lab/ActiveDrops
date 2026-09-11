#!/usr/bin/env bash
# Bounded geometry and impermeable-scalar software checks; no drop evolution.
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"
if [[ -f .project_config ]]; then
  # shellcheck disable=SC1091
  source .project_config
fi
qcc_bin="$(command -v qcc)"
if [[ "$qcc_bin" != "$repo_root/basilisk/src/qcc" ||
      ! -f basilisk/.comphy-lock ]]; then
  echo "A pinned project-local Basilisk compiler is required." >&2
  exit 1
fi
test_dir="$(mktemp -d "${TMPDIR:-/tmp}/active-drops-embed.XXXXXX")"
trap 'rm -rf -- "$test_dir"' EXIT
cp testCases/embed-contract.c "$test_dir/"
cp testCases/activity-phase-contract.c "$test_dir/"
cp testCases/activity-source-budget.c "$test_dir/"
cd "$test_dir"
for pipe in 0 1; do
  "$qcc_bin" -O2 -Wall -disable-dimensions -DTEST_PIPE="$pipe" \
    -I"$repo_root/src-local" embed-contract.c -o "embed-$pipe" -lm
  "./embed-$pipe"
done
"$qcc_bin" -O2 -Wall -disable-dimensions -I"$repo_root/src-local" \
  activity-phase-contract.c -o activity-phase-contract -lm
./activity-phase-contract
for pipe in 0 1; do
  "$qcc_bin" -O2 -Wall -disable-dimensions -DTEST_PIPE="$pipe" \
    -I"$repo_root/src-local" activity-source-budget.c \
    -o "activity-source-budget-$pipe" -lm
  "./activity-source-budget-$pipe"
done
