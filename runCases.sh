#!/bin/bash
# Compile and run the fixed-Pe demonstration (Pe = 1.6, level 8).
# For the Péclet-number scan use Script/PeScan.py; see README.md.
set -euo pipefail
cd "$(dirname "$0")"

if [[ -f .project_config ]]; then
  # shellcheck source=/dev/null
  source .project_config
fi

qcc -O2 -Wall -disable-dimensions dropMove.c -o dropMove -lm
./dropMove
