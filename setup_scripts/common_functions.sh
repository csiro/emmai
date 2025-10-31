#!/bin/bash

# Common logging function
log() {
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] $*"
}

# Common environment checking
check_environment() {
    local env_name=$1
    if $CONDA_COMMAND env list | grep -Fq "$env_name"; then
        log "Environment $env_name already exists"
        return 0
    fi
    return 1
}

# Common conda configuration setup
setup_conda_config() {
    local conda_root=$1
    mkdir -p ${conda_root}/{pkgs,envs}
    conda config --add envs_dirs $conda_root/envs
    conda config --add pkgs_dirs $conda_root/pkgs
}

# Common test data setup
setup_test_data() {
    local script_dir=$1
    local analyses_root=$2

    log "Setting up test data"
    if [ ! -d "$analyses_root" ]; then
	    mkdir -p "$analyses_root"
    fi
    cp -r "$script_dir"/../test/* "$analyses_root"
}

# Common kernel setup
setup_jupyter_kernel() {
    local env_name=$1
    local display_name=$2
    log "Setting up Jupyter kernel for $env_name"
    python -m ipykernel install --name "$env_name" --display-name "$display_name" --user
}

# Common PyTorch installation
install_pytorch() {
    local cuda_version=${1:-"cu124"}  # Default to CUDA 12.4
    log "Installing PyTorch with CUDA $cuda_version"
    pip install torch torchvision torchaudio --index-url "https://download.pytorch.org/whl/${cuda_version}"
}

# Common data download
download_model_data() {
    local unikp_root=$1
    log "Downloading model data to $unikp_root"
    mkdir -p "$unikp_root/prot_t5_xl_uniref50"
    python3 "${SCRIPT_DIR}/model_download.py" "$unikp_root"
}

# Common environment validation
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
