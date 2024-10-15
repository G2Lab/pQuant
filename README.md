# pQuant: Secure gene expression quantification from RNA-seq data

# Overview

pQuant performs secure computation on encrypted RNA-seq reads to quantify the gene expression levels for all genes by leveraging homomoprhic encryption. 
This implementation uses the OpenFHE library to ensure robust, secure operations without compromising from computational efficiency. 
Below are the detailed instructions to help you set up and run the pQuant environment.

### Analysis pipeline
The pQuant pipeline is managed using Snakemake, which handles task execution based on defined rules. The snakemake workflow for plain code is provided in the workflow/ directory. This directory includes scripts for reference creation, read concatenation, and the main pQuant steps. The pipeline replaces the previous use of Slurm scripts and optimizes execution for both local and distributed computing environments.

## Installation

### OS requirements

The package has been tested on Linux, specifically CentOS 7.9.2009. It can also run on macOS, but the CMake file may need to be revised for that. In the current version, running on Linux is recommended.

### Prerequisites and Dependencies

Please intall the following dependencies to guarantee compatibility with the pQuant setup:

 - cmake: version 3.16.3
 - gcc: version 11.2.0
 - llvm: version 9.0.0
 - clang: version 9.0.0
 - snakemake: version 7.30.1

### Installing OpenFHE

OpenFHE is a public library that supports various homomorphic encryption schemes, especially BFV scheme. Please refer to the [official documentation](https://openfhe-development.readthedocs.io/en/latest/sphinx_rsts/intro/installation/installation.html) for details. We specifically use version 1.1.2.

Follow the steps below for a successful installation on Linux:
```bash
    git clone https://github.com/openfheorg/openfhe-development.git
    cd openfhe-development
    git checkout v1.1.2
    mkdir build
    cd build
    cmake .. -DCMAKE_INSTALL_PREFIX=../../openfhe #install in local directory
    make -j
    make install
```
The installation for openFHE takes less than 10 minutes in our system.

### Building pQuant
To compile the pQuant project:
```bash
    # cd ~
    mkdir build
    cd build
    cmake ..
    make -j
```
The compilation takes less than 1 minutes in our system. Upon successful compilation, the executable named pquant is generated. Test the installation using:
```bash
    # cd ~/build
    ./pquant -t bench
```
This command will test the BFV schemes from the OpenFHE library and output runtime results for encode, encrypt, decrypt, add, and multiply operations.

## Running pQuant using Snakemake

In current version, pQuant is executed via Snakemake, which allows for management of computational resources and task dependencies.

### Configuration

Before running the pipeline, adjust the parameters in the `config/config.yaml` file:
```yaml
    K: 15
    THRES: 0.0001
    N_BATCH: 3
    exe: "build/pquant"
    OUT_DIR: "out"
    GENE_PATH: "dataset/five_genes/five_gene_reference.fa"
    READ_PATH: "dataset/five_genes/five_gene_reads.fa"

    slurm:
    step1:
        mem: "200G"
        cpus: 8
    step2:
        mem: "1G"
        cpus: 1
    step3:
        mem: "200G"
        cpus: 8
    step4:
        mem: "200G"
        cpus: 8
    step5:
        mem: "200G"
        cpus: 8
```

### Running the Pipeline

To run the entire pQuant workflow using Snakemake, execute the following commands:
```bash
    snakemake --cores <NUMBER_OF_CORES>
```

### Submitting the Snakemake Workflow to SLURM

To submit the workflow to SLURM, use the following command:
```bash
    sbatch submit_snakemake.sh
```

Here is the example of the `submit_snakemake.sh` script:
```sh
#!/bin/bash
#SBATCH --job-name=snakemake_workflow       # Job name
#SBATCH --error=slurm_out/snakemake-%j.err  # Error file for Snakemake logs
#SBATCH --output=slurm_out/snakemake-%j.out # Output file for Snakemake logs
#SBATCH --cpus-per-task=8                   # Number of CPUs for Snakemake process
#SBATCH --mem=16G                           # Memory to request for the Snakemake process

# Set job_id to the SLURM job ID
job_id=$SLURM_JOB_ID
mkdir -p out/${job_id}/slurm
# Run Snakemake with SLURM integration, using the slurm-jobscript.sh for each rule submission
snakemake -j 256 --config job_id=${job_id} \
  --cluster "sbatch --mem={resources.mem_mb}M --cpus-per-task={resources.cpus} \
    --time=12:00:00 \
    --output=out/${job_id}/slurm/{rule}-%j.out \
    --error=out/${job_id}/slurm/{rule}-%j.err" \
  --latency-wait 60 --printshellcmds --use-conda
```

### Parameters

he key parameters used in the pipeline are:
	- K: K-mer size
	- THRES: Entropy threshold for filtering
	- N_BATCH: Number of batches to divide the gene data into for parallel processing
	- GENE_PATH: Path to the reference gene file
	- READ_PATH: Path to the RNA-seq reads
	- MEM: Memory allocated for each step of the pipeline
	- OUT_DIR: Output directory where results will be saved

### Step-by-Step Execution

Snakemake automatically orchestrates the steps in the pipeline. The key steps include:
	- Step 1: Generate K-mer table (cloud)
	- Step 2: HE key generation (local)
	- Step 3: Encode and encrypt data (local)
	- Step 4: Compute matching batch (cloud)
	- Step 5: Decrypt and return gene expression vector (local)

Each step is configured with resource constraints (memory, CPUs) and logs outputs to the specified OUT_DIR.

### Results

The output will be saved in the <OUT_DIR> folder, organized by job ID and step:
	- `<OUT_DIR>/<JOB_ID>/bfv`: Stores all BFV-related files, including HE keys and ciphertexts.
	- `<OUT_DIR>/<JOB_ID>/kmer`: Stores K-mer tables and indices generated by pQuant.
	- `<OUT_DIR>/<JOB_ID>/logs`: Contains logs for each step of the Snakemake execution.
A summary of the runtime for each step will also be saved in the `time_summary.csv` file in the output directory.

This README structure reflects the new Snakemake-based workflow, simplifying job submission and resource management compared to the previous Slurm-based approach.


### Analysis pipeline
A Snakemake pipeline used to produce our results is included in the "workflow/" subdirectory. This directory also contains scripts for reference creation and read concatenation for pQuant.
