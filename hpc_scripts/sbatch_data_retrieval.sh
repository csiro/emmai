#!/bin/bash
#SBATCH --job-name=ecGEMs_dataretrieval   # Job name
#SBATCH --output=logs/SLURM-%x.%j.out     # Standard output and error log
#SBATCH --error=logs/SLURM-%x.%j.err      # Error log
#SBATCH --nodes=1                         # Run on a single Node
#SBATCH --ntasks-per-node=1               # Number of CPU cores per task
#SBATCH --mem=1GB                         # Total memory limit
#SBATCH --time=00:45:00                   # Time limit hrs:min:sec
#SBATCH --partition=io

# set up environment
. ../../env.sh
mamba activate ${ENV_ID}_cpu

export OMP_NUM_THREADS=$SLURM_NTASKS

cd ../python_scripts
python 1_data_retrieval.py
