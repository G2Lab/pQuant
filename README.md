# PQuant_HE - BFV implementation of secure pQuant algorithm

## How to run

Please follow the steps below in order.
Load required modules first
```bash
    module load cmake/3.16.3
    module load gcc/9.2.0
    module load clang
```

### Install OpenFHE

OpenFHE is a public library that supports various homomorphic encryption schemes, especially BFV scheme. Please refer [official document](https://openfhe-development.readthedocs.io/en/latest/sphinx_rsts/intro/installation/installation.html) for details. To install it on Linux, you can do as follows.
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

### Build

```bash
    mkdir build
    cd build
    cmake ..
    make -j
```
As a result, executable file is made as name `pquant`. To test, try run
```bash
    # cd ~/build
    ./pquant -t bench
```
Then the code tests BFV schemes from openFHE library, and outputs runtimes for encode/encrypt/decrypt/add/mult.

### Run
Use `*.sh` files in `job_submit/` folder. To run all jobs, run as follows
```bash
    # cd ~
    cd job_submit
    sbatch run_all.sh -k <K> -g <GENE_PATH> -r <READ_PATH> -t <THRES> -n <N_GENES> -b <N_BATCH> -m <MEM> -o <OUT_DIR>
```
`run_all.sh` runs every steps in order, batching 4th step. To change paramters, modify `run_all.sh` file from line 16. You can adjust following paramters:
- K: k
- THRES: H, threshold of entropy
- N_GENES: number of genes in dataset
- N_BATCH: number of batches you want to run. The code automatically divides genes into N_BATCH number of batches and run STEP 4 in parallel. Each job runs with {N_GENES / N_BATCH} genes
- GENE_PATH: gene(reference) path
- READ_PATH: read path
- MEM: allocated memory for each batched jobs (e.g. "200G" for 200GB memory)
- OUT_DIR: directory that all output files (include saved table, HE contents, etc) are saved

For instance, you can run with
```bash
    # cd ~
    cd job_submit
    GENE_PATH="/gpfs/commons/groups/gursoy_lab/cwalker/projects/pquant/workflow/data/reference/pquant/5k_random_protein_coding_genes.combined_exons.exons.fa"
    READ_PATH="/gpfs/commons/groups/gursoy_lab/cwalker/projects/pquant/workflow/data/test_fastqs/5k_random_protein_coding_genes.genes_only.fq"
    # note that above dataset consists of 5000 genes, so add that as an argument `-n`
    # with `-b 250`, the code automatically divides 5000 genes into 250 batches and run with 20 genes per one job in Step 4.
    sbatch run_all.sh -k 15 -t 0.00001 -g ${GENE_PATH} -r ${READ_PATH} -n 5000 -b 250 -m 200G -o ../out
```