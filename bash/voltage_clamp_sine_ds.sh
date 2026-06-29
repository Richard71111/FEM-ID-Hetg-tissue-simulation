#!/bin/bash
# =============================================================================
# Slurm array job — random sine-wave 2-cell voltage-clamp dataset
#
# Before the FIRST submission, create the log directory:
#   mkdir -p Output/save_bash_out/vc_sine_dataset
#
# Submit from the PROJECT ROOT:
#   sbatch bash/voltage_clamp_sine_ds.sh
#
# To generate more than (array_size × CASES_PER_TASK) cases, resubmit with
# a higher CASE_ID_OFFSET:
#   CASE_ID_OFFSET=10000 sbatch bash/voltage_clamp_sine_ds.sh
# =============================================================================
#SBATCH --account=PAS1622
#SBATCH --mail-type=END,FAIL
#SBATCH -J VC_sine
#SBATCH --time=24:00:00
#SBATCH --array=1-10%384        # 10 tasks × 1000 cases/task = 10 000 cases total
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --chdir=/users/PAS1622/richardsui01/FEM-ID-Hetg-tissue-simulation
#SBATCH --output=Output/save_bash_out/vc_sine_dataset/%A_%a.out
#SBATCH --error=Output/save_bash_out/vc_sine_dataset/%A_%a.err

# ----------------------------- user settings -----------------------------
CASES_PER_TASK=1000    # cases run per MATLAB process
CASE_ID_OFFSET=0       # first case_id = CASE_ID_OFFSET + 1  (override via env)
GJ_coupling="strong"
MATLAB_SCRIPT="run_sine_2CellVoltage_clamp_dataset"
# --------------------------------------------------------------------------

set -euo pipefail

# Allow CASE_ID_OFFSET to be overridden from the environment at submit time.
# Example: CASE_ID_OFFSET=10000 sbatch bash/voltage_clamp_sine_ds.sh
CASE_ID_OFFSET="${CASE_ID_OFFSET:-0}"

module load matlab/r2024b

export CASES_PER_TASK
export CASE_ID_OFFSET
export GJ_coupling
export OMP_NUM_THREADS=1

matlab -singleCompThread -batch "${MATLAB_SCRIPT}"
