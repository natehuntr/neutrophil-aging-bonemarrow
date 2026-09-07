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

if [[ -f slurm/env.sh ]]; then
  # shellcheck disable=SC1091
  source slurm/env.sh
else
  echo "WARNING: slurm/env.sh not found; falling back to the environment" >&2
fi

# Everything below runs under `set -u`, so give these definite values rather
# than referencing them bare and failing with "unbound variable".
R_MODULE="${R_MODULE:-}"
R_LIBS_USER="${R_LIBS_USER:-}"

if [[ -n "$R_MODULE" ]] && command -v module &>/dev/null; then
  module load "$R_MODULE"
fi
command -v Rscript &>/dev/null || {
  echo "ERROR: no Rscript on PATH after loading '${R_MODULE:-<none>}'." >&2
  echo "Set R_MODULE in slurm/env.sh to a module from: module -t avail | grep '^R/'" >&2
  exit 1
}

if [[ -z "$R_LIBS_USER" ]]; then
  cat >&2 <<MSG
ERROR: R_LIBS_USER is empty, so there is nowhere to install to.

Set it in slurm/env.sh, on storage the compute nodes can see and with room for
several GB -- a home directory is usually too small:

  R_LIBS_USER="\${R_LIBS_USER:-/path/to/lab/storage/R/library-\${R_MODULE//\//-}}"

Current slurm/env.sh:
$(sed -n '/^R_MODULE=/p;/^R_LIBS_USER=/p' slurm/env.sh 2>/dev/null || echo '  <file not found>')
MSG
  exit 1
fi

mkdir -p "$R_LIBS_USER" || {
  echo "ERROR: cannot create $R_LIBS_USER -- check the path and your quota." >&2
  exit 1
}
export R_LIBS_USER

echo "startup    : R_ENVIRON_USER=${R_ENVIRON_USER:-<default>} R_PROFILE_USER=${R_PROFILE_USER:-<default>}"
echo "R          : $(command -v Rscript)"
echo "             $(Rscript -e 'cat(R.version.string)')"
echo "library    : $R_LIBS_USER"
echo "cpus       : ${SLURM_CPUS_PER_TASK:-1}"
echo
echo "R's library search path (the first entry is where packages will land):"
Rscript -e 'cat(paste0("  ", .libPaths()), sep = "\n")'
echo
# A library built for another R version on this path is the single most common
# cause of an install that appears to work and then fails to load.
Rscript -e '
  bad <- grep(paste0("library/", getRversion()[, 1:2], "$"), .libPaths(),
              invert = TRUE, value = TRUE)
  bad <- grep("x86_64-pc-linux-gnu-library/[0-9]+\\.[0-9]+$", bad, value = TRUE)
  bad <- bad[basename(bad) != paste(getRversion()[, 1:2])]
  if (length(bad))
    warning("library paths built for a different R version are on the search ",
            "path: ", paste(bad, collapse = ", "), call. = FALSE)'
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
