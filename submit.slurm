#!/usr/bin/env bash
# submit.slurm -- SLURM submission TEMPLATE for the H-bond strength pipeline.
#
# This repository's production machine has no batch scheduler, so run_hb.sh is
# the primary entry point (`./run_hb.sh 4M`). Use this file only on an HPC
# cluster with SLURM. Adjust the #SBATCH lines and the module/venv block to your
# site, then submit with:
#
#   sbatch submit.slurm                 # all systems in one job
#   sbatch --export=ALL,STRIDE=100 submit.slurm     # override a default
#
# Or as a job array (one system per task) -- uncomment "#SBATCH --array" below
# and the Option-B block, then:  sbatch submit.slurm
#
#SBATCH --job-name=hb_strength
#SBATCH --output=hb_%x_%j.out
#SBATCH --error=hb_%x_%j.err
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
## #SBATCH --partition=CHANGE_ME        # <-- set your partition/queue
## #SBATCH --array=0-2                   # <-- enable for the array (Option B)

set -euo pipefail

# --- environment: adjust to your cluster ----------------------------------- #
module purge 2>/dev/null || true
# module load python/3.11                       # <-- your Python module, or:
# source "$HOME/venvs/mdanalysis/bin/activate"  # <-- your virtualenv/conda env
# The interpreter must have MDAnalysis + numpy importable.

HERE="$(cd "$(dirname "$(scontrol show job "$SLURM_JOB_ID" 2>/dev/null | \
       awk -F= '/Command=/{print $2}')" 2>/dev/null || echo .)" && pwd)"
[[ -f "$HERE/run_hb.sh" ]] || HERE="$SLURM_SUBMIT_DIR"

# --- choose which systems to run ------------------------------------------- #
# Option A (default): every system, sequentially, in this single job.
SYSTEMS_TO_RUN=(4M 4M_500 4M_1000)

# Option B: one system per array task (needs #SBATCH --array=0-2 above).
# ALL=(4M 4M_500 4M_1000)
# SYSTEMS_TO_RUN=("${ALL[$SLURM_ARRAY_TASK_ID]}")

# STRIDE (and any other config.sh default) is inherited from the environment.
bash "$HERE/run_hb.sh" "${SYSTEMS_TO_RUN[@]}"
