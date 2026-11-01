#!/bin/bash
# ======================================================================
# sim_run.sbatch
#
# Purpose  : SLURM array wrapper for R/sim_run.R. One task = one
#            (cell, block) pair = 150 fields x 30 units on one core.
# Produces : data/sim/raw/cell_<i>_block_<b>.rds, logs/<job>_<array>_<task>.out
# Reads    : R/sim_common.R, R/sim_run.R, data/sim/cells.csv
# Run      : submit through R/sim_submit.sh from the repository root on the
#            cluster login node (it creates logs/ and data/sim/raw, checks the
#            R module, and adds -p normal --requeue):
#              bash R/sim_submit.sh 1-288                     # block 1, balanced cells
#              bash R/sim_submit.sh 289-336                   # block 1, unbalanced cells
#              bash R/sim_submit.sh 1-288 --export=ALL,BLOCK=2
#              bash R/sim_submit.sh 1-$(wc -l < data/sim/tasklist_next.txt) \
#                   --export=ALL,TASKLIST=data/sim/tasklist_next.txt
#            TASKLIST: file with one "cell block" pair per line (written by
#            sim_aggregate.R); otherwise cell = array index, block = $BLOCK (1).
# Seed     : set inside sim_run.R (stream per (cell, field)). A requeued task
#            whose output already exists and matches exits without recomputing.
# ======================================================================
#SBATCH -J hx_sim
#SBATCH -p normal
#SBATCH --requeue
#SBATCH -c 1
#SBATCH --mem=2G
#SBATCH -t 02:00:00
#SBATCH -o logs/%x_%A_%a.out

set -euo pipefail

R_MODULE="R/4.5.2"
export R_LIBS_USER="$HOME/R/library"
RSCRIPT="Rscript"

module load "$R_MODULE"
cd "$SLURM_SUBMIT_DIR"
mkdir -p data/sim/raw "$R_LIBS_USER"

if [ -n "${TASKLIST:-}" ]; then
  read -r CELL BLOCK < <(sed -n "${SLURM_ARRAY_TASK_ID}p" "$TASKLIST")
else
  CELL="$SLURM_ARRAY_TASK_ID"
  BLOCK="${BLOCK:-1}"
fi
if [ -z "${CELL:-}" ] || [ -z "${BLOCK:-}" ]; then
  echo "no (cell, block) for array index ${SLURM_ARRAY_TASK_ID}" >&2; exit 2
fi

echo "task ${SLURM_ARRAY_TASK_ID}: cell ${CELL} block ${BLOCK} on $(hostname) at $(date)"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
status=0
"$RSCRIPT" R/sim_run.R "$CELL" "$BLOCK" || status=$?
if [ "$status" -ne 0 ]; then
  echo "R failed with status ${status} (cell ${CELL} block ${BLOCK}) at $(date)" >&2
  exit "$status"
fi
echo "done at $(date)"
