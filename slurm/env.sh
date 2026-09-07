# ---------------------------------------------------------------------------
# Where R comes from on this cluster. Sourced by run_pipeline.sbatch,
# submit_all.sh and install_dependencies.sbatch, so the install and the jobs
# cannot drift apart -- which is the failure mode this file exists to prevent.
#
# EDIT THESE TWO LINES for your site, then nothing else needs changing.
#
# Values already in the environment win, so a one-off override still works:
#   R_MODULE=R/4.4.0-gfbf-2023b ./slurm/submit_all.sh 7
#
# That also means re-sourcing this file after editing it does NOT pick up the
# change in a shell that already sourced it -- the old exported value wins. To
# see what the file itself says, read it in a clean environment:
#
#   env -u R_MODULE -u R_LIBS_USER bash -c 'source slurm/env.sh; echo "$R_LIBS_USER"'
#
# Submitted jobs are unaffected: they start from the submitting shell, so a
# stale value there does propagate -- unset it or open a new shell after an
# edit.
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

# ---------------------------------------------------------------------------
# Ignore personal R startup files for pipeline runs.
#
# ~/.Renviron, ~/.Rprofile and ~/.R/Makevars are read by every R session and
# override the environment the job sets up. On a cluster with several R
# versions that reliably goes wrong:
#
#   - an R_LIBS_USER pinned to another R version's library, so packages
#     install into (and load from) a tree built for a different R. The
#     symptoms are "This is R 4.3.2, package 'Matrix' needs >= 4.4",
#     "unable to load shared object .../libs/foo.so", and "failed to lock
#     directory ... for modifying";
#   - a Makevars adding a conda prefix to CPPFLAGS/LDFLAGS, so packages
#     compile against conda headers and link conda libraries while running
#     under a module R. That is what produces a dyn.load of a conda
#     libstdc++ from an R that has nothing to do with conda.
#
# This affects pipeline runs only; interactive R still reads them normally.
# Set CLEAN_R_STARTUP= (empty) to opt out.
CLEAN_R_STARTUP="${CLEAN_R_STARTUP:-1}"
if [[ -n "$CLEAN_R_STARTUP" ]]; then
  export R_ENVIRON_USER=/dev/null
  export R_PROFILE_USER=/dev/null
  export R_MAKEVARS_USER=/dev/null
  # R_LIBS would prepend other trees ahead of R_LIBS_USER.
  unset R_LIBS || true
fi
