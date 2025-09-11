# HPC Workflow Orchestration Using SLURM

This directory contains a series of SLURM batch jobs to automate a HPC workflow
of the EMMAi pipeline. The workflow consists of three jobs executed in sequence:
data retrieval, GPU-based computation, and model modification. All tasks can be
initiated using a master submission script. HPC log management and debugging
practices are covered below.

---

## Workflow Overview

The workflow involves the following SLURM batch scripts and dependencies:

1. **`sbatch_data_retrieval.sh`**
   - Take less then an hour to run for most models.
   - Retrieves sequence data from UniProt if a protein FASTA file is not
     provided, and metabolite SMILES from ChemSpider.
   - Must be executed on the **IO partition**.
   - Has inbuilt checkpointing so that you can re-run the scripts if they
     exceed their WALL time or fail for some other reason.

2. **`sbatch_uni_kp.sh`**
   - Executes the primary computation utilizing GPU resources.
   - Runs on the **GPU partition**.

3. **`sbatch_model_modifications.sh`**
   - Performs model tuning or modifications.
   - Runs on the **CPU partition**.

4. **`sbatch_submission.sh`**
   - The master script responsible for orchestrating the execution of the above
     batch scripts.
   - Configured to ensure the jobs are submitted in the correct order and
     partition constraints are respected.

You can execute the scripts via sbatch separately or you can run the shell
script which will execute them in order on the correct partitions.

---

## Prerequisites

Before running the scripts, ensure that:

- You have access to the HPC environment and necessary SLURM partitions: **IO**,
  **GPU**, and **CPU**.
- You have installed the CONDA prerequisites - best done using the install
  script [see the detailed README](../setup_scripts/README_HPC.md)
- All required modules and dependencies for your Python scripts are preloaded in
  the environment or within the batch scripts. This is most easily done by
  sourcing the env.sh file configured by the HPC installation script.

---

## File Descriptions

- **SLURM Batch Scripts**
  - `sbatch_data_retrieval.sh`: Executes the data retrieval step on the
      **IO partition**.
  - `sbatch_uni_kp.sh`: Executes GPU-bound processing on the **GPU partition**.
  - `sbatch_model_modifications.sh`: Runs CPU-based model modifications on the
      **CPU partition**.
  
- **Master Submission Script**
  - `sbatch_submission.sh`: Facilitates sequential submission of the above
      steps in the correct order and partitions.

- **Log Directory**
  - All SLURM logs are saved to a designated directory (`logs/`) to
      assist in debugging.

---

## Step-by-Step Instructions

### 1. Master Script: `sbatch_submission.sh`

The master script `sbatch_submission.sh` is used to orchestrate the workflow. It
ensures that:

- The first job (`sbatch_data_retrieval.sh`) runs on the IO partition and
  completes successfully before the next step is initiated.
- The second job (`sbatch_uni_kp.sh`) runs on the GPU partition after the first
  job completes. This will be sent to the GPU partition without any local
  environment variables from the CPU cluster, except what is explicitly set
  in the sbatch command. (This may be important to consider when debugging)
- The third job (`sbatch_model_modifications.sh`) runs on the CPU partition
  after the second job completes.

**Execution:**
Navigate to the directory containing the scripts and run:

```bash
./sbatch_submission.sh
```

Note: each script will source the env.sh file from the root directory of your
installation.

Note: there is an example multi-sample script `sbatch_multi_submission.sh` in
case you have multiple analyses. It constrains the analysis so that only one IO
process is run at a time in an attempt to not overload remote APIs.
