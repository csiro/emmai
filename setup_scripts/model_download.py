#!/usr/bin/python3

import sys
from huggingface_hub import hf_hub_download
from typing import List

def download_model_files(download_path: str, repo_id: str, files: List[str]) -> None:
    """
    Download model files from Hugging Face Hub.

    Args:
        download_path (str): Path where files will be cached
        repo_id (str): Hugging Face repository ID
        files (List[str]): List of files to download
    """
    for file in files:
        hf_hub_download(
            repo_id=repo_id,
            filename=file,
            local_dir=download_path
        )

# Ensure command line argument is provided
if len(sys.argv) != 2:
    print("Usage: script.py DATA_PATH")
    sys.exit(1)

# Files to download
data_path = sys.argv[1]
prot_t5_files = [
    "pytorch_model.bin",
    "tokenizer_config.json",
    "config.json",
    "spiece.model",
    "special_tokens_map.json"
]
unikp_files = ["UniKP for kcat.pkl"]

# Download files
download_model_files(f"{data_path}/prot_t5_xl_uniref50", "Rostlab/prot_t5_xl_uniref50", prot_t5_files)
download_model_files(data_path, "HanselYu/UniKP", unikp_files)
