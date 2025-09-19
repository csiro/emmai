#!/bin/bash
#SBATCH --job-name=emmai_model_mod       # Job name
#SBATCH --output=logs/SLURM-%x.%j.out     # Standard output and error log
#SBATCH --error=logs/SLURM-%x.%j.err      # Error log
#SBATCH --nodes=1                         # Run on a single Node
#SBATCH --ntasks-per-node=1               # Number of CPU cores per task
#SBATCH --mem=1GB                         # Total memory limit
#SBATCH --time=00:05:00                   # Time limit hrs:min:sec
#SBATCH --partition=defq

# set up environment
. ../../env_emmai_cpu.sh
mamba activate emmai_cpu

export OMP_NUM_THREADS=$SLURM_NTASKS

cd ../python_scripts
python 3_model_modification.py
echo "Finished 3_model_modification"
python 4_patching_models.py
echo "Finished 4_patching_models"
python 5_protein_pool_calibration.py
echo "Finished 5_patching_models"
