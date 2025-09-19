#!/bin/bash
#SBATCH --job-name=emmai_uni_kp          # Job name
#SBATCH --output=logs/SLURM-%x.%j.out     # Standard output and error log
#SBATCH --error=logs/SLURM-%x.%j.err      # Error log
#SBATCH --nodes=1                         # Run on a single Node
#SBATCH --ntasks-per-node=4               # Number of CPU cores per task
#SBATCH --mem=2GB                         # Total memory limit
#SBATCH --time=00:10:00                   # Time limit hrs:min:sec
#SBATCH --gpus-per-node=2
#SBATCH --partition=gpu

# set up environment
. ../../env_emmai_gpu.sh
mamba activate emmai_gpu

export OMP_NUM_THREADS=$SLURM_NTASKS

cd ../python_scripts
python 2_uni_kp_prot.py
