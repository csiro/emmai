# EMMAi

## Description

The `emmai` repository contains a collection of Jupyter notebooks, Python scripts, and Slurm batch files designed to facilitate the EMMAi computational workflow.

## Contents

- **Jupyter Notebooks**: Interactive notebooks for data analysis and visualization.
- **Python Scripts**: Standalone scripts for various computational tasks. These for the most part duplicate the code in the Jupyter notebooks.
- **Slurm Batch Files**: Scripts to manage and execute jobs on Slurm-based high-performance computing clusters. These execute the Python scripts after activating a Conda environment.

## Requirements

To use the contents of this repository, you will need the following:
- Mamba (or Conda if you prefer)

## Installation

1. Clone the repository:

   ```bash
   git clone --recurse-submodules git@github.com:csiro-internal/emmai.git
   ```

2. There are two installation scripts provided in the setup_scripts directory depending on where you intend to run the pipeline.
- HPC systems with GPU and CPU clusters - [see the detailed README](setup_scripts/README_HPC.md)
- Standalone systems - [see the detailed README](setup_scripts/README_standalone.md)

## Executing the pipeline

After following the installation guides you can choose to execute the
pipeline in one of three ways, or a combination of the three - if you understand how each method works.

1. By running the Jupyter notebooks sequentially starting with:
    ```
    notebooks/1_data_retrieval.ipynb
    ```

2. By executing the python scripts sequentially starting with:
    ```
    python_scripts/1_data_retrieval.py
    ```

3. Or by running a combination of the batch scripts in hpc_scripts. See [See the HPC README](hpc_scripts/README.md) for details.
