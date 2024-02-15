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
    mkdir build
    cd build
    cmake .. -DCMAKE_INSTALL_PREFIX=../openfhe #install in local directory
    make
    make install
```
After that, test library as
```bash
    make testall
```


### Build

```bash
    mkdir build
    cd build
    cmake ..
```


### Run
Use `*.sh` files in `job_submit/` folder. To run all jobs, run as follows
```bash
    # cd ~
    cd job_submit
    sbatch run_all.sh
```
`run_all.sh` runs every steps in order, batching 4th step. To change paramters, modify `run_all.sh` file from line 16. You can adjust following paramters:
- K: k
- THRES: H, threshold of entropy
- N_GENES: number of genes in dataset
- N_BATCH: number of batches you want to run. The code automatically divides genes into N_BATCH number of batches and run STEP 4 in parallel. Each job runs with {N_GENES / N_BATCH} genes
- GENE_PATH: gene(reference) path
- READ_PATH: read path
- MEM: allocated memory for each batched jobs (e.g. "200G" for 200GB memory)


Our algorithm consists of 5 steps. For instance, run step 1 as
```bash
    # cd build
     ./pquant -t STEP1 -k 5 --thres 0.01 -g {gene_path} -r {read_path}
```

