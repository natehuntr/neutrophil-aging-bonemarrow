#!/bin/bash
# ---------------------------------------------------------------------------
# Submit the pipeline as one job per step, chained by SLURM dependencies.
#
#   ./slurm/submit_all.sh              # submit everything
#   ./slurm/submit_all.sh 5 7 9        # only these steps (dependencies still
#                                      # apply between the ones you name)
#   DRY_RUN=1 ./slurm/submit_all.sh    # print the sbatch commands, submit none
#
# Why not one big job: the steps want very different resources, a failure only
# costs the step that failed, and the graph allows genuine parallelism.
#
#   1 -> 2 -> 3 -+-> 4
#                +-> 6
#                +-> 8
#                +-> 5 -+-> 7
#                       +-> 9
#
# Steps 4, 6 and 8 start together once 3 lands; 7 and 9 once 5 does. Each job
# runs with `afterok`, so a failure stops that branch instead of feeding a
# later step a half-written object.
# ---------------------------------------------------------------------------

set -euo pipefail

cd "$(dirname "$0")/.."
[[ -f config/config.yml ]] || { echo "run this from the project root" >&2; exit 1; }
mkdir -p logs

SBATCH_SCRIPT=slurm/run_pipeline.sbatch
REQUESTED=("$@")

# Per-step resources. These are STARTING POINTS, not measurements: check
# `seff <jobid>` after the first successful run and tighten them. Step 7 is the
# outlier because it refits every model glm_de.n_perm times.
resources_for() {
  case "$1" in
    1) echo "--time=04:00:00 --mem=64G  --cpus-per-task=4" ;;
    2) echo "--time=12:00:00 --mem=64G  --cpus-per-task=4" ;;
    3) echo "--time=08:00:00 --mem=96G  --cpus-per-task=4" ;;
    4) echo "--time=04:00:00 --mem=48G  --cpus-per-task=4" ;;
    5) echo "--time=24:00:00 --mem=64G  --cpus-per-task=8" ;;
    6) echo "--time=12:00:00 --mem=48G  --cpus-per-task=4" ;;
    7) echo "--time=48:00:00 --mem=64G  --cpus-per-task=8" ;;
    8) echo "--time=06:00:00 --mem=32G  --cpus-per-task=4" ;;
    9) echo "--time=02:00:00 --mem=32G  --cpus-per-task=4" ;;
  esac
}

# Which already-submitted steps each step must wait for.
depends_on() {
  case "$1" in
    1) echo "" ;;
    2) echo "1" ;;
    3) echo "2" ;;
    4|5|6|8) echo "3" ;;
    7|9) echo "5" ;;
  esac
}

wanted() {
  [[ ${#REQUESTED[@]} -eq 0 ]] && return 0
  local step
  for step in "${REQUESTED[@]}"; do [[ "$step" == "$1" ]] && return 0; done
  return 1
}

declare -A JOB_ID=()

for step in 1 2 3 4 5 6 7 8 9; do
  wanted "$step" || continue

  dep_args=""
  dep_list=""
  for dep in $(depends_on "$step"); do
    if [[ -n "${JOB_ID[$dep]:-}" ]]; then
      dep_list="${dep_list}${dep_list:+:}${JOB_ID[$dep]}"
    fi
  done
  [[ -n "$dep_list" ]] && dep_args="--dependency=afterok:${dep_list}"

  # shellcheck disable=SC2046
  cmd=(sbatch --parsable --job-name="bm-step${step}" $(resources_for "$step")
       ${dep_args:+$dep_args} --export=ALL,STEPS="$step" "$SBATCH_SCRIPT")

  if [[ -n "${DRY_RUN:-}" ]]; then
    printf '%s\n' "${cmd[*]}"
    JOB_ID[$step]="<step${step}>"
  else
    JOB_ID[$step]=$("${cmd[@]}")
    printf 'step %s -> job %s%s\n' "$step" "${JOB_ID[$step]}" \
      "${dep_list:+  (after ${dep_list})}"
  fi
done

if [[ -z "${DRY_RUN:-}" ]]; then
  cat <<'NOTE'

Submitted. Useful afterwards:
  squeue -u "$USER" -o '%.10i %.12j %.9T %.10M %R'   # what is queued/running
  tail -f logs/bm-step3-<jobid>.out                  # follow a step
  seff <jobid>                                       # what it actually used
  scancel -u "$USER" --name=bm-step7                 # cancel one step

A step that fails leaves its dependents in DependencyNeverSatisfied; fix the
cause, then resubmit that step and the ones after it.
NOTE
fi
