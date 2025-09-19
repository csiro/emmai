# EMMAi HPC Environment Setup

This directory contains scripts to set up the `emmai` HPC environment,
including installing necessary dependencies, cloning required repositories,
and setting up the directory structure.

## HPC Installation Scripts and Files

- `install.sh`: Automates the setup of the environment, including creating directories, cloning repositories, and downloading data.
- `hpc_cpu_requirements.yml`: Environment configuration of the Conda CPU environment.
- `hpc_gpu_requirements.yml`: Environment configuration of the Conda GPU environment.

## Installation

### Prerequisites

#### Software
Make sure you have the following installed on your system:

- Conda or Mamba
- Jupyter (for notebooks)
- Slurm workload manager (for batch files)

#### Disk space
- Ensure that you have sufficient permissions and at least 20GB disk space for the project directories and model files.

### Setting Up the Environment

To set up the environment, follow these steps:

#### Step 1: Update the PROJECT_ROOT path

Edit the `install.sh` file to set your project root path and your CPU and GPU hosts. Note if you only have a CPU or GPU cluster you can set the same hostname for both of these variables. For example:

```env
PROJECT_ROOT=/home/<your_username>/emmai
CPU_HOST=petrichor-login
GPU_HOST=virga-login
```

The above assumes your PROJECT_ROOT directory is in your home directory, but it is probably more appropriate to place this on a high performance scratch file system on your cluster if this is available.

#### Step 2: Run the `install.sh` Script

On each cluster (if you have a CPU and GPU cluster), navigate to the `setup_scripts` directory if not already there and run the `install.sh`
script and specify the `hpc` mode:

```bash
cd setup_scripts
bash install.sh hpc
```

This script will:

1. Set up the directory structure.
2. Create and activate the CONDA environment
3. Install additional packages using pip and create a Jupyter kernel.
4. Download required model data files.

#### Step 3: Verify the Setup

After the script completes, you can activate the environment by sourcing the `env_<env_name>.sh` file in the `PROJECT_ROOT` directory. The below assumes you placed your `PROJECT_ROOT` in your home directory.

```bash
source /home/<your_username>/emmai/env_<env_name>.sh
```
