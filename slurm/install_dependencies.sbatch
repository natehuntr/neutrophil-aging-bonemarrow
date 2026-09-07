#!/bin/bash
# ---------------------------------------------------------------------------
# Install the pipeline's R packages once, into the library named in
# slurm/env.sh. Run this before the first pipeline submission.
#
#   sbatch slurm/install_dependencies.sbatch
#
# A batch job rather than an interactive session because most of these compile
# from source and it takes a while -- submit it and come back.
# ---------------------------------------------------------------------------

#SBATCH --job-name=bm-install
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --nodes=1
#SBATCH --ntasks=1

set -euo pipefail

cd "${PROJECT_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
[[ -f config/config.yml ]] || { echo "ERROR: run from the project root" >&2; exit 1; }
mkdir -p logs

# shellcheck disable=SC1091
[[ -f slurm/env.sh ]] && source slurm/env.sh

if [[ -n "${R_MODULE:-}" ]] && command -v module &>/dev/null; then
  module load "$R_MODULE"
fi
command -v Rscript &>/dev/null || { echo "ERROR: no Rscript on PATH" >&2; exit 1; }

mkdir -p "$R_LIBS_USER"

echo "R          : $(command -v Rscript)"
echo "             $(Rscript -e 'cat(R.version.string)')"
echo "library    : $R_LIBS_USER"
echo "cpus       : ${SLURM_CPUS_PER_TASK:-1}"
echo

Rscript scripts/00_install_dependencies.R
status=$?

echo
echo "install finished with status $status"
if [[ $status -eq 0 ]]; then
  echo "Verifying the pipeline can load what it needs:"
  Rscript -e '
    source("R/packages.R")
    ok <- tryCatch({ load_packages(); TRUE }, error = function(e) { message(conditionMessage(e)); FALSE })
    quit(status = as.integer(!ok))'
  status=$?
fi
exit $status
