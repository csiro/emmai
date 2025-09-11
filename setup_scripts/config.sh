#!/bin/bash
export PROJECT_ROOT=~/ecgems
export CPU_HOST=petrichor-login
export GPU_HOST=virga-login
export ANALYSIS_ROOT=${PROJECT_ROOT}/analyses
export CONDA_ROOT=${PROJECT_ROOT}/conda
export ENV_ID=ecGEMs
export PYTHON_VERSION=3.9
export SCIKIT_LEARN_VERSION=0.24.2
# git clone --recurse-submodules
export UNIKP=${PROJECT_ROOT}/ecgems_pipeline/unikp
if [ -z "${PYTHONPATH:-}" ]; then
    export PYTHONPATH=$UNIKP
else
    export PYTHONPATH=$UNIKP:$PYTHONPATH
fi
# Use mamba if installed
if command -v mamba &> /dev/null; then
    export CONDA_COMMAND=mamba
else
    export CONDA_COMMAND=conda
fi