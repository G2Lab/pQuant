# PQuant: Secure gene expression quantification from RNA-seq data

# Overview

pQuant performs secure computation on encrypted RNA-seq reads to quantify the gene expression levels for all genes by leveraging homomoprhic encryption. 
This implementation uses the OpenFHE library to ensure robust, secure operations without compromising from computational efficiency. 
Below are the detailed instructions to help you set up and run the pQuant environment.

## Installation

### Prerequisites and Dependencies

Please intall the following dependencies to guarantee compatibility with the pQuant setup:

- cmake: version 3.16.3 or later
- gcc: version 9.2.0 or later
- llvm: version 9.0.0 or later
- clang: version 9.0.0 or later

### Installing OpenFHE

OpenFHE is a public library that supports various homomorphic encryption schemes, especially BFV scheme. Please refer to the [official document](https://openfhe-development.readthedocs.io/en/latest/sphinx_rsts/intro/installation/installation.html) for details. 
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

### Building pQuant
To compile the pQuant project:
```bash
    mkdir build
    cd build
    cmake ..
    make -j
```
Upon successful compilation, the executable named pquant is generated. Test the installation using:
```bash
    # cd ~/build
    ./pquant -t bench
```
This command will test the BFV schemes from the OpenFHE library and output runtime results for encode, encrypt, decrypt, add, and multiply operations.

## Execution via Slurm

pQuant is optimized for batch processing using Slurm. Scripts are located in the `job_submit/` folder.

### Running All Jobs
Navigate to the job submission directory and execute the jobs with:
```bash
    cd ~/job_submit
    sbatch run_all.sh -k <K> -g <GENE_PATH> -r <READ_PATH> -t <THRES> -n <N_GENES> -b <N_BATCH> -m <MEM> -o <OUT_DIR>
```

The parameters we use are as listed:
- K: k, size of k-mer
- THRES: H, threshold of entropy
- N_GENES: number of genes in dataset
- N_BATCH: number of batches you want to run. The code automatically divides genes into N_BATCH number of batches and run STEP 4 in parallel. Each job runs with {N_GENES / N_BATCH} genes
- GENE_PATH: gene(reference) path (.fa)
- READ_PATH: read path (.fa, .fq, or .fastq)
- MEM: allocated memory for each batched jobs (e.g. "200G" for 200GB memory)
- OUT_DIR: directory that all output files (include saved table, HE contents, etc) are saved

### Results

The results are stored in `<OUT_DIR>` folder. The slurm script creates a folder of name `<SLURM_JOB_ID>` and three other folders `bfv`, `kmer`, and `slurm_out`
- `<OUT_DIR>/<SLURM_JOB_ID>/bfv`: all bfv-related stuffs are stored, including context, HE keys, and ciphertexts
- `<OUT_DIR>/<SLURM_JOB_ID>/kmer`: all results from pQuant algorithm are stored. In our implementation, kmer table and the index of kmer are stored.
- `<OUT_DIR>/<SLURM_JOB_ID>/slurm_out`: all logs are stored. 
