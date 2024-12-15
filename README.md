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

### Install with Docker

You can use Docker to simplify the setup and run pQuant without manually installing dependencies. Follow these steps:

1. Install docker: install docker for your operating system.
2.Build the Docker Image

The Dockerfile is located in the Docker/ directory. Use the following command to build the Docker image:
```bash
    docker build -t pquant:latest -f Docker/Dockerfile .
```
3. Run pQuant Using Docker

After building the image, you can run pQuant inside a Docker container.
Start the container interactively to test the pQuant commands:
```bash
    docker run -it pquant:latest /bin/bash
```

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
snakemake -j 256 --config job_id=${job_id} OUT_DIR=${OUT_DIR} \
  --cluster "sbatch --mem={resources.mem_mb}M --cpus-per-task={resources.cpus} --job-name={rule} \
    --time=48:00:00 \
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

## Summary Statistics

To address user requests and improve the usability of pQuant, we have introduced additional functionalities that generate summary statistics at key stages of the pipeline. These enhancements provide insights into input data quality and final output, ensuring data reliability and usability while maintaining privacy guarantees. All these steps are required to take less than a minute.

### Pre-encoding Statistics

`Note`: seqtk/1.2 must be available and loaded in your environment to enable this feature.

Before the encoding stage, pQuant generates a set of summary statistics from the input RNA-seq reads:
- Total number of k-mers: The number of k-mers identified in the input dataset.
- Total number of reads: The count of all reads in the dataset.
- Number of reads meeting the quality threshold: Using the seqtk library, this functionality filters reads based on a user-defined quality threshold.

These statistics allow users to verify input data quality and assess its suitability for downstream processing.

This part is commanded in snakemake `pre_analysis` rule. To run with snakemake and toy dataset, please follow following:
```bash
    # Ensure seqtk/1.2 is loaded in your environment
    module load seqtk/1.2

    # Run the pQuant algorithm. The following command runs the end-to-end process on the test dataset,
    # but for pre-analysis, you only need to execute step 1.
    snakemake -j 8 --configfile config/config_toy.yaml

    # The job_id is automatically set in the format YYMMDD-HHMMSS, as specified in the config_toy.yaml file.
    # Check your output directory for the generated job_id before running the following command with your <JOB_ID>.
    # You can manage the quality control threshold using the <qc> parameter.

    snakemake -j 8 --configfile config/config_toy.yaml --config job_id=<JOB_ID> QC=<qc> --printshellcmds --rerun-incomplete pre_analysis
```
For instance, the output from toy dataset is as following in `data_summary/pre_analysis.txt'.
```text
    ...    
    === read reads ===
    read file is .fa format
    readFastaFile duration = 58 ms
    Total number of kmers: 98140
    Total number of reads: 5
    Number of reads in out.fa with quality > 5: 80453

```

### Post-decryption Statistics

After the decryption stage, pQuant provides metrics on the gene expression results:
- Number of genes with TPM values above a user-defined threshold: Users can specify a TPM (Transcripts Per Million) threshold, and the functionality calculates the number of genes meeting this criterion.

These metrics allow users to evaluate the completeness and reliability of the quantification results.

This step is implemented in the post_analysis rule of the Snakemake pipeline. To execute it, follow the steps below:
```bash
    # Run the pQuant algorithm.
    snakemake -j 8 --configfile config/config_toy.yaml

    # The job_id is automatically set in the format YYMMDD-HHMMSS, as specified in the config_toy.yaml file.
    # Check your output directory for the generated job_id before running the following command with your <JOB_ID>.
    # You can manage the TPM threshold parameter using the <thres> variable.

    snakemake -j 8 --configfile config/config_toy.yaml --config job_id=<JOB_ID> TPM_THRES=<thres> --printshellcmds --rerun-incomplete post_analysis
```
For instance, the output from toy dataset is as following in `data_summary/post_analysis.txt'.
```text
    ...    
    Gene: ENSG00000131174, TPM: 382830
    Gene: ENSG00000119707, TPM: 129092
    Gene: ENSG00000105971, TPM: 9478.97
    Gene: ENSG00000089289, TPM: 446375
    Gene: ENSG00000085982, TPM: 32224.1
    Number of genes with TPM > 0.5: 5
```