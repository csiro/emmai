# ECGEMS Standalone Environment Setup

This directory contains scripts to set up the `ecgems` standalone environment,
including installing necessary dependencies, cloning required repositories,
and setting up the directory structure.

## Repository Contents

- `standalone_environment.yml`: Environment configuration, including
  dependencies and channels for the Conda environment.
- `standalone_install.sh`: Automates the setup of the environment,
  including creating directories, cloning repositories, and downloading data.
- `config.sh`: Contains variables used to by the `standalone_install.sh` script

## Installation

### Prerequisites

Make sure you have the following installed on your system:

- Conda or Mamba
- Git

### Cloning the repository

```bash
git clone --recurse-submodules https://github.com/csiro-internal/emmai.git
```

### Setting Up the Environment

To set up the environment, follow these steps:

#### Step 1: Update the PROJECT_ROOT path

Navigate to the `setup_scripts` directory if not already there and edit the `config.sh` file to set your project root path.
For example:

```env
PROJECT_ROOT=/home/<your_username>/ecgems
```

#### Step 2: Run the `standalone_install.sh` Script

Run the `standalone_install.sh` script to set up the environment:

```bash
bash standalone_install.sh
```

This script will:

1. Set up the directory structure.
2. Create and activate the Conda environment using the
   `standalone_environment.yml` file.
3. Install additional packages using pip and create a Jupyter kernel.
4. Download required model data files.

#### Step 3: Verify the Setup

After the script completes, verify that the environment has been activated.

## Notes

- Ensure that you have sufficient permissions and at least 20GB disk space for
  the project directories and model files.
