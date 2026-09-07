#!/bin/bash
# ---------------------------------------------------------------------------
# Run the pipeline (or some of its steps) as one SLURM job.
#
#   sbatch slurm/run_pipeline.sbatch                    # every step, 1-9
#   sbatch --export=ALL,STEPS="5" slurm/run_pipeline.sbatch
#   sbatch --export=ALL,STEPS="6 8 9" slurm/run_pipeline.sbatch
#
# Resources here cover the heaviest step in a full run. To give one step its
# own allocation, override at submit time rather than editing this file:
#
#   sbatch --time=48:00:00 --mem=64G --export=ALL,STEPS="7" \
#          slurm/run_pipeline.sbatch
#
# Submit from the project root: SLURM_SUBMIT_DIR is where the scripts expect
# to start, and the job checks for config/config.yml before doing anything.
# ---------------------------------------------------------------------------

#SBATCH --job-name=bm-aging
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --nodes=1
#SBATCH --ntasks=1

# --- EDIT THESE FOR YOUR CLUSTER -------------------------------------------
# Left commented so the job submits against your site defaults as-is. Uncomment
# and fill in whichever your cluster requires.
##SBATCH --partition=general
##SBATCH --account=your_account
##SBATCH --qos=normal
##SBATCH --mail-type=END,FAIL
##SBATCH --mail-user=you@example.edu

# How to get R. Pick ONE; setting both is what breaks a conda setup, because a
# module's R ends up first on PATH while the libraries still come from conda.
#
#   R_MODULE   name of an environment module to load, e.g. "R/4.3.1"
#   CONDA_ENV  name or path of a conda environment to activate
#
# Both default to empty, in which case whatever R is on PATH at submit time is
# used. Jobs are submitted with --export=ALL, so a conda environment activated
# before sbatch carries through on its own and needs nothing set here.
R_MODULE="${R_MODULE:-}"
CONDA_ENV="${CONDA_ENV:-}"

# Set to 1 to run `module purge` first. Off by default: it strips compiler and
# library modules that an R build may depend on.
MODULE_PURGE="${MODULE_PURGE:-}"

# Optional personal library. Leave unset to use R's default.
# export R_LIBS_USER="$HOME/R/library"
# ---------------------------------------------------------------------------

set -euo pipefail

# Find the project root: PROJECT_DIR if given, otherwise walk up from wherever
# the job was submitted. Walking up means submitting from slurm/ or scripts/
# works as well as from the root. $0 is no help here -- SLURM runs a copy of
# this script from a spool directory, so it says nothing about where the
# project is.
find_project_root() {
  local dir
  dir="$(cd "${1:-.}" 2>/dev/null && pwd)" || return 1
  while [[ "$dir" != "/" ]]; do
    [[ -f "$dir/config/config.yml" ]] && { printf '%s\n' "$dir"; return 0; }
    dir="$(dirname "$dir")"
  done
  return 1
}

START_DIR="${PROJECT_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
if ! PROJECT_ROOT="$(find_project_root "$START_DIR")"; then
  echo "ERROR: no config/config.yml in $START_DIR or any directory above it." >&2
  echo "" >&2
  echo "Submit from inside the project, or name it explicitly:" >&2
  echo "  cd /path/to/the/project && sbatch slurm/run_pipeline.sbatch" >&2
  echo "  sbatch --export=ALL,PROJECT_DIR=/path/to/the/project slurm/run_pipeline.sbatch" >&2
  echo "" >&2
  echo "To locate it:" >&2
  echo "  find \"$START_DIR\" -maxdepth 3 -path '*/config/config.yml'" >&2
  exit 1
fi
cd "$PROJECT_ROOT"
mkdir -p logs

if [[ -n "$R_MODULE" && -n "$CONDA_ENV" ]]; then
  echo "ERROR: set R_MODULE or CONDA_ENV, not both." >&2
  echo "Loading a module on top of a conda environment gives you the module's" >&2
  echo "R running against conda's libraries, which fails at dyn.load." >&2
  exit 1
fi

if [[ -n "$R_MODULE" ]]; then
  command -v module &>/dev/null || { echo "ERROR: R_MODULE is set but this node has no 'module' command." >&2; exit 1; }
  [[ -n "$MODULE_PURGE" ]] && { module purge || true; }
  module load "$R_MODULE"
fi

if [[ -n "$CONDA_ENV" ]]; then
  # `conda activate` needs the shell hook; it is not available in a
  # non-interactive shell by default.
  CONDA_BASE="$(conda info --base 2>/dev/null || true)"
  if [[ -z "$CONDA_BASE" || ! -f "$CONDA_BASE/etc/profile.d/conda.sh" ]]; then
    echo "ERROR: CONDA_ENV is set but conda was not found on this node." >&2
    echo "Either load your conda module first, or activate the environment" >&2
    echo "before sbatch and leave CONDA_ENV unset -- --export=ALL carries it." >&2
    exit 1
  fi
  # shellcheck disable=SC1091
  source "$CONDA_BASE/etc/profile.d/conda.sh"
  conda activate "$CONDA_ENV"
fi

if ! command -v Rscript &>/dev/null; then
  echo "ERROR: Rscript is not on PATH." >&2
  echo "Set R_MODULE, or CONDA_ENV, or activate your environment before sbatch." >&2
  exit 1
fi

# R's linear algebra will otherwise grab every core on the node and thrash
# against the allocation. Keep BLAS single-threaded and let the pipeline's own
# parallelism (compute.cores, used by monocle3) do the work.
CPUS="${SLURM_CPUS_PER_TASK:-1}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export R_DATATABLE_NUM_THREADS="$CPUS"

STEPS="${STEPS:-}"

echo "=========================================================="
echo "job          : ${SLURM_JOB_ID:-<interactive>} on $(hostname)"
echo "project root : $PWD"
echo "steps        : ${STEPS:-all (1-9)}"
echo "cpus / mem   : $CPUS / ${SLURM_MEM_PER_NODE:-?} MB"
echo "R            : $(command -v Rscript)"
echo "               $(Rscript -e 'cat(R.version.string)' 2>/dev/null)"
[[ -n "${CONDA_PREFIX:-}" ]] && echo "conda env    : $CONDA_PREFIX"
echo "R_LIBS_USER  : ${R_LIBS_USER:-<unset>}"
echo "started      : $(date)"
echo "=========================================================="

# compute.cores in config/config.yml is what monocle3 is told to use. Warn if
# the allocation cannot back it.
CONFIG_CORES=$(grep -E '^\s+cores:' config/config.yml | head -1 | grep -oE '[0-9]+' || echo 1)
if (( CONFIG_CORES > CPUS )); then
  echo "WARNING: compute.cores is $CONFIG_CORES but this job has $CPUS CPUs."
  echo "         Lower it in config/config.yml or raise --cpus-per-task."
fi

# `set -e` would abort here before the footer runs, losing both the timing and
# the exit code, so the failure is caught explicitly instead.
set +e
# shellcheck disable=SC2086
Rscript scripts/run_all.R $STEPS
status=$?
set -e

echo "=========================================================="
echo "finished     : $(date) (exit $status)"
echo "=========================================================="
exit $status
