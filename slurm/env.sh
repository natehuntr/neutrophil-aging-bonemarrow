# ---------------------------------------------------------------------------
# Where R comes from on this cluster. Sourced by run_pipeline.sbatch,
# submit_all.sh and install_dependencies.sbatch, so the install and the jobs
# cannot drift apart -- which is the failure mode this file exists to prevent.
#
# EDIT THESE TWO LINES for your site, then nothing else needs changing.
#
# Values already in the environment win, so a one-off override still works:
#   R_MODULE=R/4.4.0-gfbf-2023b ./slurm/submit_all.sh 7
# ---------------------------------------------------------------------------

# An environment module providing R. Leave empty if R is already on PATH (e.g.
# a conda environment activated before submitting).
R_MODULE="${R_MODULE:-R/4.3.2-gfbf-2023a}"

# Where the pipeline's own packages live. Two requirements:
#   - on a filesystem the compute nodes can see (NOT a node-local /tmp, and on
#     some clusters not $HOME either);
#   - versioned by R release, so switching modules does not mix incompatible
#     builds in one directory.
#
# $HOME is the default but is often the wrong choice: home quotas are small and
# this library runs to several GB. Prefer group or lab storage, e.g.
#   R_LIBS_USER="${R_LIBS_USER:-/nemo/lab/<lab>/home/users/<user>/R/library-${R_MODULE//\//-}}"
R_LIBS_USER="${R_LIBS_USER:-$HOME/R/library-${R_MODULE//\//-}}"

export R_MODULE R_LIBS_USER
