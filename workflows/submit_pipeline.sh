#!/bin/bash
# Submit full transpoly pipeline with job dependencies

set -euo pipefail

WORKFLOWS_DIR="$(dirname "$0")"

# Submit equilibration job
EQ_JOB=$(sbatch "$WORKFLOWS_DIR/run_equilibration.slurm" | awk '{print $4}')
echo "Submitted equilibration job: $EQ_JOB"

# Submit production job, dependent on equilibration
PROD_JOB=$(sbatch --dependency=afterok:"$EQ_JOB" "$WORKFLOWS_DIR/run_production.slurm" | awk '{print $4}')
echo "Submitted production job: $PROD_JOB (depends on $EQ_JOB)"

echo "Pipeline submitted. Check status with: squeue -j $EQ_JOB,$PROD_JOB"
