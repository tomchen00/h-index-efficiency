#!/bin/bash
# ======================================================================
# sim_driver.sh
#
# Purpose  : Runs the remaining rounds of the simulation stopping rule on the cluster login node:
#            waits for the current array job, aggregates (as a batch job), and submits the next
#            round from data/sim/tasklist_next.txt until no cell needs another block or the
#            block cap (4 blocks = 600 fields) is reached; then writes a final aggregation.
# Run      : from the repository root on the login node --
#              nohup bash R/sim_driver.sh <current_jobid> <current_round> > logs/driver.log 2>&1 &
# ======================================================================
set -uo pipefail
JOB="$1"; ROUND="$2"; MAXROUND=4
wait_job() {
  while [ "$(squeue -h -r -j "$1" 2>/dev/null | wc -l)" -gt 0 ]; do sleep 60; done
  sacct -j "$1" -n -X --format=State%20 | sort | uniq -c
}
aggregate() {
  local a
  a=$(sbatch --parsable -p quicktest -t 00:15:00 --mem=32G -c 1 -J hx_agg -o logs/agg_%j.out \
        --wrap="module load R/4.5.2; Rscript R/sim_aggregate.R $*")
  echo "$(date) aggregation job $a"
  wait_job "$a"
  tail -3 "logs/agg_${a}.out"
}
while true; do
  echo "$(date) waiting for round $ROUND job $JOB"
  wait_job "$JOB"
  aggregate
  N=$(wc -l < data/sim/tasklist_next.txt 2>/dev/null || echo 0)
  if [ "$N" -eq 0 ] || [ "$ROUND" -ge "$MAXROUND" ]; then
    echo "$(date) stopping: round $ROUND, $N cells still above the MCSE target"
    break
  fi
  ROUND=$((ROUND + 1))
  cp data/sim/tasklist_next.txt "data/sim/tasklist_round${ROUND}.txt"
  JOB=$(bash R/sim_submit.sh "1-${N}" --export=ALL,TASKLIST="data/sim/tasklist_round${ROUND}.txt" | grep -o 'job [0-9]*' | awk '{print $2}')
  echo "$(date) round $ROUND submitted: job $JOB, $N tasks" | tee -a submitted.log
done
echo "$(date) DRIVER DONE"
