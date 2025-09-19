# EMMAi Standalone Environment Setup

This directory contains scripts to set up the `emmai` standalone environment, including installing necessary dependencies, cloning required repositories, and setting up the directory structure.

## Repository Contents

- `install.sh`: Automates the setup of the environment, including creating directories, cloning repositories, and downloading data.
- `standalone_environment.yml`: Environment configuration, including dependencies and channels for the Conda environment.

## Installation

### Prerequisites

#### Software
Make sure you have the following installed on your system:

- Conda or Mamba
- Jupyter (for notebooks)

#### Disk space
- Ensure that you have sufficient permissions and at least 20GB disk space for the project directories and model files.
### Setting Up the Environment

To set up the environment, follow these steps:

#### Step 1: Update the PROJECT_ROOT path
Navigate to the `setup_scripts` directory if not already there and edit the `install.sh` file to set your project root path.
For example:

```env
PROJECT_ROOT=/home/<your_username>/emmai
```

#### Step 2: Run the `install.sh` Script
Run the `install.sh` script with the `standalone` mode to set up the environment:

```bash
bash install.sh standalone
```

This script will:

1. Set up the directory structure.
2. Create and activate the Conda environment using the `standalone_environment.yml` file.
3. Install additional packages using pip and create a Jupyter kernel.
4. Download required model data files.

#### Step 3: Verify the Setup

After the script completes, you can activate the environment with the command provided by the script, for example:
```bash
Activate environment with: conda activate /home/kubeflow/emmai_project/conda
```
