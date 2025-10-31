#!/bin/bash
set -e
SCRIPT_DIR=$(dirname "$(realpath "$0")")
SMILES_DB=/scratch3/projects/ecgems/DATABASES/SMILES_reference_DB.csv
source "${SCRIPT_DIR}/config.sh"
source "${SCRIPT_DIR}/common_functions.sh"

setup_hpc_gpu_environment() {
    local env_id=$1
    local python_version=$2
    local scikit_learn_version=$3

    log "Setting up HPC GPU environment"
    conda config --add channels conda-forge
    conda config --add channels bioconda
    mamba create -y -n "${env_id}_gpu" python=$python_version cmake swig scikit-learn=$scikit_learn_version
    mamba activate "${env_id}_gpu"
    install_pytorch "cu124"
    mamba install -y --file "${SCRIPT_DIR}/hpc_gpu_requirements.txt"
}

setup_hpc_cpu_environment() {
    local env_id=$1
    local python_version=$2

    log "Setting up HPC CPU environment"
    conda config --add channels conda-forge
    conda config --add channels bioconda
    mamba create -y -n "${env_id}_cpu" python=$python_version
    mamba activate "${env_id}_cpu"
    mamba install -y --file "${SCRIPT_DIR}/hpc_cpu_requirements.txt"
}

# Setup based on host
HOST=$(hostname)
if [ "$HOST" == "$CPU_HOST" ]; then
    if [ -f ${SMILES_DB} ]; then
        cp ${SMILES_DB} ${SCRIPT_DIR}/../python_scripts
    else
        echo "No SMILES reference Database found."
        echo "Copy SMILES_reference_DB.csv to the Python scripts directory"
    fi
    module purge
    module load miniforge3
    MAMBA_EXE=$(which mamba)
    . ${MINIFORGE3_HOME}/etc/profile.d/conda.sh
    . ${MINIFORGE3_HOME}/etc/profile.d/mamba.sh
    setup_conda_config "$CONDA_ROOT"
    validate_environment || exit 1
    if ! check_environment "${ENV_ID}_cpu"; then
        setup_test_data "$SCRIPT_DIR" "$ANALYSES_ROOT"
        setup_hpc_cpu_environment "$ENV_ID" "$PYTHON_VERSION"
        mamba activate ${ENV_ID}_cpu
        setup_jupyter_kernel "$ENV_ID_cpu" "$ENV_ID_cpu"
    fi
fi
if [ "$HOST" == "$GPU_HOST" ]; then
    module purge
    module load miniforge3
    MAMBA_EXE=$(which mamba)
    . ${MINIFORGE3_HOME}/etc/profile.d/conda.sh
    . ${MINIFORGE3_HOME}/etc/profile.d/mamba.sh
    setup_conda_config "$CONDA_ROOT"
    validate_environment || exit 1
    if ! check_environment "${ENV_ID}_gpu"; then
        setup_hpc_gpu_environment "$ENV_ID" "$PYTHON_VERSION" "$SCIKIT_LEARN_VERSION"
        mamba activate ${ENV_ID}_gpu
        setup_jupyter_kernel "$ENV_ID_gpu" "$ENV_ID_gpu"
        download_model_data "$UNIKP"
    fi
fi

cat <<EOF >${PROJECT_ROOT}/env.sh
#!/bin/bash
source "${SCRIPT_DIR}/config.sh"
# Change the INPUTS variable to the directory where your inputs.yml file diamond
# model file and your optional protein fasta file are located.
# Do not set this if submitting multiple analyses to the cluster
export INPUTS=\${ANALYSES_ROOT}/PAO1

module purge
module load miniforge3
MAMBA_EXE=\$(which mamba)
. \${MINIFORGE3_HOME}/etc/profile.d/conda.sh
. \${MINIFORGE3_HOME}/etc/profile.d/mamba.sh

HOST=\$(hostname)
if [ "\$HOST" == "\$CPU_HOST" ]; then
    mamba activate \${ENV_ID}_cpu
elif [ "\$HOST" == "\$GPU_HOST" ]; then
    mamba activate \${ENV_ID}_gpu
fi
EOF

log "Activate environment with source ${PROJECT_ROOT}/env.sh"
