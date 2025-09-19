#!/usr/bin/env bash
set -eu
trap 'echo "Script interrupted"; exit 1' INT TERM

PROJECT_ROOT=/scratch3/wat440/emmai_project
CPU_HOST=petrichor-i3
GPU_HOST=virga-i2
CPU_CONDA_MODULE="miniconda3/23.3.1"
GPU_CONDA_MODULE="miniconda3/23.5.2"
ANALYSIS_ROOT="${PROJECT_ROOT}/analyses"
CONDA_ROOT="${PROJECT_ROOT}/conda"
UNIKP="../unikp"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODE="$1"

# Logging
log() {
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] $*"
}

# Check if Conda environment exists
check_environment() {
    local env_path="$1"
    if [ -d "$env_path" ] && [ -d "$env_path/conda-meta" ]; then
        log "Environment at $env_path already exists"
        return 0
    fi
    return 1
}

# Create Conda environment
setup_environment() {
    local env_id="$1"
    local env_yaml="$2"
    local env_path="$3"
    local conda_command

    # Use mamba if installed
    if command -v mamba &> /dev/null; then
        conda_command=mamba
    else
        conda_command=conda
    fi

    log "Setting up '$env_id' environment"
    if ! $conda_command env create -f "$env_yaml" --prefix "$env_path"; then
        log "Failed to create environment from $env_yaml"
        exit 1
    fi

    log "Installing PyTorch with CUDA cu124"
    "${ENV_PATH}/bin/pip" install torch torchvision torchaudio --index-url "https://download.pytorch.org/whl/cu124"

    log "Setting up Jupyter kernel for $env_id"
    "${ENV_PATH}/bin/python" -m ipykernel install --name "$env_id" --display-name "$env_id" --user
}
# Load HPC modules
init_modules() {
    local conda_module="$1"
    module purge
    module load "$conda_module"
    . "${MINICONDA3_HOME}/etc/profile.d/conda.sh"
    . "${MINICONDA3_HOME}/etc/profile.d/mamba.sh"
}

# Validate required commands
validate_environment() {
    local required_commands=("conda" "pip" "python")
    for cmd in "${required_commands[@]}"; do
        if ! command -v "$cmd" >/dev/null 2>&1; then
            log "Error: Required command '$cmd' not found"
            return 1
        fi
    done
    return 0
}

# Usage help
usage() {
    echo "Usage: $0 <mode>"
    echo "  mode: 'standalone' or 'hpc'"
    echo "Example:"
    echo "  $0 standalone"
    echo "  $0 hpc"
    exit 1
}

# Entry point
if [[ $# -ne 1 || "$MODE" == "--help" || "$MODE" == "-h" ]]; then
    usage
fi

if [ -z "${PYTHONPATH:-}" ]; then
    export PYTHONPATH="$UNIKP"
else
    export PYTHONPATH="$UNIKP:$PYTHONPATH"
fi

validate_environment || exit 1
mkdir -p "$PROJECT_ROOT" "$ANALYSIS_ROOT"

# Mode-specific setup
if [ "$MODE" == "standalone" ]; then
    if command -v conda &>/dev/null; then
        . "$(conda info --base)/etc/profile.d/conda.sh"
    else
        log "Conda is not installed or not in PATH"
        exit 1
    fi
    ENV_NAME="emmai"
    ENV_PATH="$CONDA_ROOT"

    if ! check_environment "$ENV_PATH"; then
        setup_environment "$ENV_NAME" "$SCRIPT_DIR/standalone_environment.yml" "$ENV_PATH"
    else
        log "Environment $ENV_PATH already exists."
    fi

    END_MESSAGE="Activate environment with: conda activate $ENV_PATH"

elif [ "$MODE" == "hpc" ]; then
    HOST=$(hostname)
    case "$HOST" in
        "$CPU_HOST")
            CONDA_MODULE="$CPU_CONDA_MODULE"
            ENV_NAME="emmai_cpu"
            ENV_FILE="hpc_cpu_environment.yml"
            ;;
        "$GPU_HOST")
            CONDA_MODULE="$GPU_CONDA_MODULE"
            ENV_NAME="emmai_gpu"
            ENV_FILE="hpc_gpu_environment.yml"
            ;;
        *)
            log "Error: Unknown HPC host '$HOST'. Expected '$CPU_HOST' or '$GPU_HOST'."
            exit 1
            ;;
    esac

    init_modules "$CONDA_MODULE"
    ENV_PATH="${CONDA_ROOT}/${ENV_NAME}"

    if ! check_environment "$ENV_PATH"; then
        setup_environment "$ENV_NAME" "$SCRIPT_DIR/$ENV_FILE" "$ENV_PATH"
    else
        log "Environment $ENV_PATH already exists."
    fi

    cat <<EOF >"${PROJECT_ROOT}/env_${ENV_NAME}.sh"
#!/usr/bin/env bash
# Optional: set INPUTS to the directory where your inputs.yml file,
# diamond model file and your optional protein fasta file are located if submitting a single analysis
# export INPUTS=\${ANALYSIS_ROOT}/analyses/PA01

module purge
module load $CONDA_MODULE
. \${MINICONDA3_HOME}/etc/profile.d/conda.sh
. \${MINICONDA3_HOME}/etc/profile.d/mamba.sh

export UNIKP="$(realpath "${SCRIPT_DIR}/${UNIKP}")"

mamba activate $ENV_PATH
EOF

    END_MESSAGE="Activate environment with: . ${PROJECT_ROOT}/env_${ENV_NAME}.sh"

else
    log "Error: Unknown mode '$MODE'"
    usage
fi

log "Downloading model data to $UNIKP"
mkdir -p "$UNIKP/prot_t5_xl_uniref50"
"${ENV_PATH}/bin/python" "${SCRIPT_DIR}/model_download.py" "$UNIKP"

log "Setting up test data"
cp -R "${SCRIPT_DIR}/../test/"* "$ANALYSIS_ROOT"

log "$ENV_NAME setup complete."
log "$END_MESSAGE"
