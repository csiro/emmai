# ecGEMs_pipeline

## Description

The `ecGEMs_pipeline` repository contains a collection of IPython notebooks,
Python scripts, and SLURM batch files designed to facilitate the ecGEMs
computational workflow.

## Contents

- **IPython Notebooks**: Interactive notebooks for data analysis and visualization.
- **Python Scripts**: Standalone scripts for various computational tasks. These
for the most part duplicate the code in the IPython notebooks.
- **SLURM Batch Files**: Scripts to manage and execute jobs on SLURM-based
high-performance computing clusters. These execute the Python scripts after
activating a CONDA environment.

## Requirements

To use the contents of this repository, you will need the following:

- Python 3.9
- Jupyter Notebook
- SLURM workload manager (for batch files)
- MAMBA (or CONDA if you prefer)

## Installation

1. Clone the repository:

   ```bash
   git clone --recurse-submodules git@github.com:csiro-internal/emmai-unikp.git
   ```

2. There are two installation scripts provided in the setup_scripts directory. 1
One for HPC systems with GPU and CPU clusters - [see the detailed README](setup_scripts/README_HPC.md) and
one for standalone systems - [see the detailed README](setup_scripts/README_standalone.md)

## Executing the pipeline

1. After following the installation guidelines you can choose to execute the
pipeline in 1 of three ways, or a combination of the three - if you know what
you are doing.

    i. By running the Jupyter notebooks sequentially starting with:
    notebooks/1_data_retrieval.ipynb

    ii. By executing the python scripts sequentially starting with:
    python_scripts/1_data_retrieval.py

    iii. Or by running a combination of the batch scripts in hpc_scripts. See
    [See the HPC README](hpc_scripts/README.md) for details.
