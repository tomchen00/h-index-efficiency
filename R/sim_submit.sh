#!/bin/bash
# ======================================================================
# sim_submit.sh
#
# Purpose  : Submits the simulation array from the cluster login node:
#            creates logs/ and data/sim/raw (SLURM opens the log file before
#            the job script runs), checks that the R module exists and loads,
#            checks the inputs, then submits R/sim_run.sbatch on
#            normal with --requeue and prints the job id.
# Produces : a SLURM array job; its id on stdout
# Reads    : R/sim_run.sbatch, R/sim_run.R, R/sim_common.R, data/sim/cells.csv
# Run      : from the repository root --
#              bash R/sim_submit.sh <array-range> [extra sbatch arguments]
#            e.g. bash R/sim_submit.sh 1-288
#                 bash R/sim_submit.sh 1-288 --export=ALL,BLOCK=2
# Seed     : n/a (seeding is inside sim_run.R)
# ======================================================================
set -euo pipefail

if [ $# -lt 1 ]; then
  echo "usage: bash R/sim_submit.sh <array-range> [extra sbatch arguments]" >&2
  exit 2
fi
RANGE="$1"; shift
R_MODULE="R/4.5.2"

mkdir -p logs data/sim/raw
for f in R/sim_run.sbatch R/sim_run.R R/sim_common.R data/sim/cells.csv; do
  [ -f "$f" ] || { echo "missing $f (run from the repository root)" >&2; exit 1; }
done

if ! command -v module >/dev/null 2>&1; then
  for init in /etc/profile.d/lmod.sh /etc/profile.d/modules.sh /usr/share/lmod/lmod/init/bash; do
    [ -f "$init" ] && source "$init" && break
  done
fi
command -v module >/dev/null 2>&1 || { echo "environment modules not available" >&2; exit 1; }
if ! module avail "$R_MODULE" 2>&1 | grep -q "$R_MODULE"; then
  echo "module $R_MODULE not found (module avail)" >&2; exit 1
fi
( module load "$R_MODULE" && command -v Rscript >/dev/null ) || { echo "module load $R_MODULE failed" >&2; exit 1; }

JOB=$(sbatch --parsable -p normal --requeue --array="$RANGE" "$@" R/sim_run.sbatch)
echo "submitted job ${JOB} (array ${RANGE})"
