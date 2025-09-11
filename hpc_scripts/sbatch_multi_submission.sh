#!/bin/bash
#
# Get the directory of the script
SCRIPT_DIR="$(dirname "$(realpath "$0")")"
cd $SCRIPT_DIR

module purge
module load SC slurm

mkdir -p logs

SPECIES1=PAO1
SPECIES2=iML1515
SPECIES3=iJO1366
declare -A JOBS

# CHANGE THIS
ANALYSES_ROOT=/scratch3/ben324/ECGEMS/analyses

export INPUTS=${ANALYSES_ROOT}/${SPECIES1}
# Assume we are on the CPU cluster in this first set of jobs each job is dependent on the previous one finishing
jobs["${SPECIES1}_1"]=$(sbatch --job-name=${SPECIES1}-IO sbatch_data_retrieval.sh | awk '{ print $4 }')
jobs["${SPECIES1}_2"]=$(sbatch --job-name=${SPECIES1}-GPU --export=NONE,INPUTS --dependency=afterok:${jobs[${SPECIES1}_1]} sbatch_uni_kp.sh | awk '{ print $4 }')
sbatch --job-name=${SPECIES1}-CPU --dependency=afterok:${jobs[${SPECIES1}_2]} sbatch_model_modifications.sh

# CHANGE THIS
export INPUTS=${ANALYSES_ROOT}/${SPECIES2}
# Assume we are on the CPU cluster and wait for first IO job
jobs["${SPECIES2}_1"]=$(sbatch --job-name=${SPECIES2}-IO --dependency=afterok:${jobs[${SPECIES1}_1]} sbatch_data_retrieval.sh | awk '{ print $4 }')
jobs["${SPECIES2}_2"]=$(sbatch --job-name=${SPECIES2}-GPU --export=NONE,INPUTS --dependency=afterok:${jobs[${SPECIES2}_1]} sbatch_uni_kp.sh | awk '{ print $4 }')
sbatch --job-name=${SPECIES2}-CPU --dependency=afterok:${jobs[${SPECIES2}_2]} sbatch_model_modifications.sh

# CHANGE THIS
export INPUTS=${ANALYSES_ROOT}/${SPECIES3}
# Assume we are on the CPU cluster and wait for second IO job
jobs["${SPECIES3}_1"]=$(sbatch --job-name=${SPECIES3}-IO --dependency=afterok:${jobs[${SPECIES2}_1]} sbatch_data_retrieval.sh | awk '{ print $4 }')
jobs["${SPECIES3}_2"]=$(sbatch --job-name=${SPECIES3}-GPU --export=NONE,INPUTS --dependency=afterok:${jobs[${SPECIES3}_1]} sbatch_uni_kp.sh | awk '{ print $4 }')
sbatch --job-name=${SPECIES3}-CPU --dependency=afterok:${jobs[${SPECIES3}_2]} sbatch_model_modifications.sh
