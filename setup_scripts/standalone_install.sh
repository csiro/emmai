#!/bin/bash
set -eu

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"
source "${SCRIPT_DIR}/common_functions.sh"

setup_standalone_environment() {
    local env_id=$1
    local env_yaml=$2

    log "Setting up '$env_id' environment"
    $CONDA_COMMAND env create -f "$env_yaml" --prefix "$CONDA_ROOT"
    conda activate "$CONDA_ROOT"
    install_pytorch "cu124"
    pip install -r "$UNIKP/requirements.txt"
    setup_jupyter_kernel "$env_id" "$env_id"
}

# Check if 'conda' is initialized
if command -v conda &> /dev/null; then
    CONDA_BASE=$(conda info --base)
    source "$CONDA_BASE/etc/profile.d/conda.sh"
else
    echo "Conda is not installed or not in PATH"
    exit 1
fi

setup_conda_config "$CONDA_ROOT"

# Validate environment
validate_environment || exit 1

# Create analyses directory
mkdir -p "$PROJECT_ROOT" "$ANALYSES_ROOT"

# Setup environment
if ! check_environment "$CONDA_ROOT"; then
    setup_standalone_environment "$ENV_ID" "$SCRIPT_DIR/standalone_environment.yml"
else
    log "Environment $CONDA_ROOT already exists."
    conda activate "$CONDA_ROOT"
fi

# Setup directories and data
setup_test_data "$SCRIPT_DIR" "$PROJECT_ROOT"
download_model_data "$UNIKP"

log "$ENV_ID setup complete."
log "Activate environment with \`conda activate $CONDA_ROOT\`"
