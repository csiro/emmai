# EMMAi HPC Environment Setup

This directory contains scripts to set up the `emmai` HPC environment,
including installing necessary dependencies, cloning required repositories,
and setting up the directory structure.

## HPC Installation Scripts and Files

- `hpc_install.sh`: Automates the setup of the environment, including creating
  directories, cloning repositories, and downloading data.
- `hpc_cpu_requirements.txt`: Environment configuration of the CONDA CPU
  environment.
- `hpc_gpu_requirements.txt`: Environment configuration of the CONDA GPU
  environment.
- `config.sh`: Environment variables that will be used by the Installation
  shell scripts.

## Installation

### Prerequisites

Make sure you have the following installed on your system and in your path:

- CONDA or Mamba: module load miniforge3

### Cloning the repository

```bash
git clone --recurse-submodules https://github.com/csiro-internal/emmai.git
```

### Setting Up the Environment

To set up the environment, follow these steps:

#### Step 1: Update the PROJECT_ROOT path

Edit the `config.sh` file to set your project root path and your CPU and GPU
hosts. Note if you only have a GPU cluster you can set the same hostname for
both of these variables. For example:

```env
export PROJECT_ROOT=/home/<your_username>/emmai
export CPU_HOST=petrichor-login
export GPU_HOST=virga-login
```

The above assumes your PROJECT_ROOT directory is in your home directory, but it
is probably more appropriate to place this on a high performance scratch file
system on your cluster if this is available.

#### Step 2: Run the `hpc_install.sh` Script

On each cluster (if you have a CPU and GPU cluster), navigate to the
`setup_scripts` directory if not already there and run the `hpc_install.sh`
script:

```bash
cd setup_scripts
bash hpc_install.sh
```

This script will:

1. Set up the directory structure.
2. Create and activate the CONDA environment
3. Install additional packages using pip and create a Jupyter kernel.
4. Download required model data files.

#### Step 3: Verify the Setup

After the script completes, you can activate the environment by sourcing the
env.sh file in the PROJECT_ROOT directory. The below assumes you placed your
PROJECT_ROOT in your home directory.

```bash
source /home/<your_username>/emmai/env.sh
```

## Notes

- Ensure that you have sufficient permissions and at least 20GB disk space for
  the project directories and model files.
