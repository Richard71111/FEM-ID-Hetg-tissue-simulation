#!/bin/bash
#
# Submit the 2-cell voltage-clamp dataset as Slurm array jobs.
# Cluster limit: MaxArraySize = 1001  ->  max array index = 1000.
# We therefore cap the array at 1000 tasks and scale total output with
# CASES_PER_TASK (cases per process) and multiple batches (CASE_ID_OFFSET).
#
# Total cases produced = TOTAL_CASES (auto-split into batches of <=1000 tasks).
#
set -euo pipefail

# ----------------------------- user settings -----------------------------
TOTAL_CASES=10000         # how many cases you want in total
CASES_PER_TASK=1000       # cases run per MATLAB process (memory: ~tens of MB each)
ARRAY_LIMIT=384          # max array tasks running at once (the %N throttle)
WALLTIME="24:00:00"      # wall time PER task (one task = CASES_PER_TASK cases)
MAX_ARRAY_INDEX=1000     # = MaxArraySize - 1 ; do not exceed on this cluster
GJ_coupling="strong"

MATLAB_SCRIPT="run_2CellVoltage_clamp_dataset"
# --------------------------------------------------------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

RUN_ID="$(date +%Y%m%d_%H%M%S)"
LOG_DIR="${PROJECT_DIR}/Output/save_bash_out/vc_dataset/FEM1/GJ_${GJ_coupling}/run_${RUN_ID}"
mkdir -p "${LOG_DIR}"

cases_per_batch=$(( MAX_ARRAY_INDEX * CASES_PER_TASK ))
num_batches=$(( (TOTAL_CASES + cases_per_batch - 1) / cases_per_batch ))

echo "RUN_ID            : ${RUN_ID}"
echo "TOTAL_CASES       : ${TOTAL_CASES}"
echo "CASES_PER_TASK    : ${CASES_PER_TASK}"
echo "cases per batch   : ${cases_per_batch}  (= ${MAX_ARRAY_INDEX} tasks x ${CASES_PER_TASK})"
echo "number of batches : ${num_batches}"
echo "Logs folder       : ${LOG_DIR}"
echo

done_cases=0
for (( b=0; b<num_batches; b++ )); do
    offset=${done_cases}
    remaining=$(( TOTAL_CASES - done_cases ))

    # tasks needed for the remaining cases, capped at the array limit
    tasks_needed=$(( (remaining + CASES_PER_TASK - 1) / CASES_PER_TASK ))
    array_end=$(( tasks_needed < MAX_ARRAY_INDEX ? tasks_needed : MAX_ARRAY_INDEX ))

    echo "Batch ${b}: array=1-${array_end}%${ARRAY_LIMIT}, CASE_ID_OFFSET=${offset}"

    sbatch <<EOT
#!/bin/bash
#SBATCH --account=PAS1622
#SBATCH --mail-type=END,FAIL
#SBATCH -J VC_ds_b${b}
#SBATCH --time=${WALLTIME}
#SBATCH --array=1-${array_end}%${ARRAY_LIMIT}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=${LOG_DIR}/vc_b${b}_%A_%a.out
#SBATCH --error=${LOG_DIR}/vc_b${b}_%A_%a.err

set -euo pipefail

module load matlab/r2024b

export CASES_PER_TASK=${CASES_PER_TASK}
export CASE_ID_OFFSET=${offset}
export RUN_ID=${RUN_ID}
export GJ_coupling=${GJ_coupling}
export OMP_NUM_THREADS=1

cd "${PROJECT_DIR}"

matlab -singleCompThread -batch "${MATLAB_SCRIPT}"
EOT

    done_cases=$(( done_cases + array_end * CASES_PER_TASK ))
done

echo
echo "Submitted ${num_batches} batch(es), covering case_id 1 .. ${done_cases}."