#!/bin/bash
#
FAST_FILE_SYSTEM=$SCRATCH3DIR
export PROJECT_ROOT=$FAST_FILE_SYSTEM/EMMAI
export CPU_HOST=petrichor-login
export GPU_HOST=virga-login
export ANALYSIS_ROOT=${PROJECT_ROOT}/analyses
export CONDA_ROOT=${PROJECT_ROOT}/conda
export ENV_ID=emmai
export PYTHON_VERSION=3.9
export SCIKIT_LEARN_VERSION=0.24.2
export UNIKP=${PROJECT_ROOT}/emmai/unikp
if [ -z "${PYTHONPATH:-}" ]; then
    export PYTHONPATH=$UNIKP
else
    export PYTHONPATH=$UNIKP:$PYTHONPATH
fi
# Use mamba if installed
if command -v mamba &>/dev/null; then
    export CONDA_COMMAND=mamba
else
    export CONDA_COMMAND=conda
fi
