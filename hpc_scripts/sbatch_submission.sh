#!/bin/bash
#
# Get the directory of the script
SCRIPT_DIR="$(dirname "$(realpath "$0")")"
cd $SCRIPT_DIR

module purge
module load SC slurm

mkdir -p logs

# Assume we are on the CPU cluster
job_id_1=$(sbatch sbatch_data_retrieval.sh | awk '{ print $4 }')
job_id_2=$(sbatch --export=NONE --dependency=afterok:$job_id_1 sbatch_uni_kp.sh | awk '{ print $4 }')
sbatch --dependency=afterok:$job_id_2 sbatch_model_modifications.sh
