#!/bin/bash
# runSimulation.sh
#
# Run a single ActiveDrops case from the repository root. The script creates
# simulationCases/c<CaseNo>/, copies the parameter file and the source file,
# compiles the case against the project-local Basilisk and runs it with
# case.params as the only argument to the binary.
#
# Usage:
#   bash runSimulation.sh [params_file] [--exec source.c] [--threads N]

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

usage() {
  cat <<'EOF'
Usage: bash runSimulation.sh [params_file] [OPTIONS]

Arguments:
  params_file    Parameter file path (default: default.params)

Options:
  --exec FILE    C source in simulationCases/ (default: dropMove.c)
  --threads N    OpenMP thread count; N=1 runs serial (default: 1)
  -h, --help     Show this help message
EOF
}

get_param_value() {
  local key="$1"
  local file="$2"
  awk -F '=' -v key="$key" '
    /^[[:space:]]*#/ { next }
    {
      k = $1
      gsub(/^[[:space:]]+|[[:space:]]+$/, "", k)
      if (k == key) {
        v = $2
        sub(/[[:space:]]*#.*/, "", v)
        gsub(/^[[:space:]]+|[[:space:]]+$/, "", v)
        print v
        exit
      }
    }
  ' "$file"
}

EXEC_CODE="dropMove.c"
PARAM_FILE="default.params"
PARAM_FILE_SET=0
OMP_THREADS=1

while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help)
      usage
      exit 0
      ;;
    --exec)
      if [[ -z "${2:-}" ]]; then
        echo "ERROR: --exec requires a file name." >&2
        usage
        exit 1
      fi
      EXEC_CODE="$2"
      shift 2
      ;;
    --exec=*)
      EXEC_CODE="${1#*=}"
      shift
      ;;
    --threads)
      if [[ -z "${2:-}" ]]; then
        echo "ERROR: --threads requires a positive integer value." >&2
        usage
        exit 1
      fi
      OMP_THREADS="$2"
      shift 2
      ;;
    --threads=*)
      OMP_THREADS="${1#*=}"
      shift
      ;;
    --)
      shift
      break
      ;;
    -*)
      echo "ERROR: Unknown option: $1" >&2
      usage
      exit 1
      ;;
    *)
      if [[ $PARAM_FILE_SET -eq 0 ]]; then
        PARAM_FILE="$1"
        PARAM_FILE_SET=1
        shift
      else
        echo "ERROR: Unexpected argument: $1" >&2
        usage
        exit 1
      fi
      ;;
  esac
done

if [[ $# -gt 0 ]]; then
  echo "ERROR: Unexpected trailing arguments: $*" >&2
  usage
  exit 1
fi

if [[ ! "$OMP_THREADS" =~ ^[1-9][0-9]*$ ]]; then
  echo "ERROR: --threads must be a positive integer, got: $OMP_THREADS" >&2
  exit 1
fi

USE_OPENMP=0
if [[ "$OMP_THREADS" -gt 1 ]]; then
  USE_OPENMP=1
fi

if [[ "$EXEC_CODE" != *.c ]]; then
  EXEC_CODE="${EXEC_CODE}.c"
fi

if [[ "$EXEC_CODE" == */* ]]; then
  echo "ERROR: --exec must be a file name inside simulationCases/, got: $EXEC_CODE" >&2
  exit 1
fi

if [[ ! "$PARAM_FILE" = /* ]]; then
  PARAM_FILE="${SCRIPT_DIR}/${PARAM_FILE}"
fi

if [[ -f "${SCRIPT_DIR}/.project_config" ]]; then
  # shellcheck disable=SC1091
  source "${SCRIPT_DIR}/.project_config"
fi

if ! command -v qcc >/dev/null 2>&1; then
  echo "ERROR: qcc not found in PATH." >&2
  echo "Hint: install the project-local Basilisk (see README.md) so that .project_config exists." >&2
  exit 1
fi

if [[ "$EXEC_CODE" == "dropMove-embed-pipe.c" ||
      "$EXEC_CODE" == "dropMove-embed-channel.c" ]]; then
  if [[ "$(command -v qcc)" != "${SCRIPT_DIR}/basilisk/src/qcc" ||
        ! -f "${SCRIPT_DIR}/basilisk/.comphy-lock" ]] ||
     ! grep -qx 'ref=v2026-08-30' "${SCRIPT_DIR}/basilisk/.comphy-lock"; then
    echo "ERROR: Embedded cases require project-local Basilisk v2026-08-30." >&2
    exit 1
  fi
fi

if [[ ! -f "$PARAM_FILE" ]]; then
  echo "ERROR: Parameter file not found: $PARAM_FILE" >&2
  exit 1
fi

if grep -Eq '^[[:space:]]*Oh[[:space:]]*=' "$PARAM_FILE"; then
  echo "ERROR: Parameter 'Oh' is retired; supply Re and Ca instead." >&2
  exit 1
fi

SRC_FILE_ORIG="${SCRIPT_DIR}/simulationCases/${EXEC_CODE}"
if [[ ! -f "$SRC_FILE_ORIG" ]]; then
  echo "ERROR: Source file not found: $SRC_FILE_ORIG" >&2
  exit 1
fi

CASE_NO="$(get_param_value "CaseNo" "$PARAM_FILE")"
if [[ -z "$CASE_NO" ]]; then
  echo "ERROR: CaseNo not found in parameter file: $PARAM_FILE" >&2
  exit 1
fi

if [[ ! "$CASE_NO" =~ ^[0-9]+$ ]]; then
  echo "ERROR: CaseNo must be numeric, got: $CASE_NO" >&2
  exit 1
fi

if [[ "$CASE_NO" -lt 1000 ]]; then
  echo "ERROR: CaseNo must be >= 1000 for consistent sorting, got: $CASE_NO" >&2
  exit 1
fi

CASE_TAG="c${CASE_NO}"
CASE_DIR="${SCRIPT_DIR}/simulationCases/${CASE_TAG}"
EXECUTABLE_NAME="${EXEC_CODE%.c}"

echo "========================================="
echo "ActiveDrops - Single Case Runner"
echo "========================================="
echo "Source file: ${EXEC_CODE}"
echo "Parameter file: ${PARAM_FILE}"
echo "CaseNo: ${CASE_NO}"
echo "Case directory: simulationCases/${CASE_TAG}"
echo "qcc: $(command -v qcc)"
if [[ $USE_OPENMP -eq 1 ]]; then
  echo "Run mode: OpenMP (threads=${OMP_THREADS})"
else
  echo "Run mode: Serial"
fi
echo "========================================="
echo ""

mkdir -p "$CASE_DIR"
cp "$PARAM_FILE" "$CASE_DIR/case.params"
cp "$SRC_FILE_ORIG" "$CASE_DIR/$EXEC_CODE"

cd "$CASE_DIR"

echo "Compiling ${EXEC_CODE} ..."
QCC_FLAGS=(-I../../src-local -O2 -Wall -disable-dimensions)
if [[ $USE_OPENMP -eq 1 ]]; then
  QCC_FLAGS+=(-fopenmp)
fi
if ! qcc "${QCC_FLAGS[@]}" "$EXEC_CODE" -o "$EXECUTABLE_NAME" -lm; then
  if [[ $USE_OPENMP -eq 1 ]]; then
    echo "ERROR: OpenMP build failed. Re-run with --threads 1 for serial mode." >&2
  fi
  exit 1
fi
echo "Compilation successful: $EXECUTABLE_NAME"
echo ""

if [[ $USE_OPENMP -eq 1 ]]; then
  echo "Running: OMP_NUM_THREADS=${OMP_THREADS} ./${EXECUTABLE_NAME} case.params"
  if OMP_NUM_THREADS="$OMP_THREADS" ./"$EXECUTABLE_NAME" case.params; then
    EXIT_CODE=0
  else
    EXIT_CODE=$?
  fi
else
  echo "Running (serial): ./${EXECUTABLE_NAME} case.params"
  if ./"$EXECUTABLE_NAME" case.params; then
    EXIT_CODE=0
  else
    EXIT_CODE=$?
  fi
fi

echo ""
if [[ $EXIT_CODE -eq 0 ]]; then
  echo "Simulation completed."
  echo "Output location: simulationCases/${CASE_TAG}/ (snapshots in intermediate/, diagnostics in log.dat)"
else
  echo "Simulation exited with code: $EXIT_CODE"
fi

exit "$EXIT_CODE"
