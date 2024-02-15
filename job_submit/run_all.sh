#!/bin/bash
#SBATCH --job-name=pq_all    # Job name
#SBATCH --partition=pe2   # Specify the partition
#SBATCH --nodes=1                 # Number of nodes
#SBATCH --ntasks-per-node=1       # Number of tasks (MPI processes) per node
#SBATCH --cpus-per-task=1         # Number of CPU cores per task
#SBATCH --mem=40G                  # Memory per node (adjust as needed)
#SBATCH --mail-user=shong@nygenome.org   # Specify your email address

# Navigate to the directory containing your executable
cd /gpfs/commons/groups/gursoy_lab/shong/Github/pQuant_enc/job_submit

echo "SLURM_JOB_ID: ${SLURM_JOB_ID}"
# Define paths and arguments

# DATA="five"
# K=10
# THRES=2
# N_GENES=5
# N_BATCH=5
# MEM="5G"

GENE_PATH=/gpfs/commons/groups/gursoy_lab/cwalker/projects/pquant/workflow/data/reference/pquant/5k_random_protein_coding_genes.combined_exons.exons.fa
READ_PATH="/gpfs/commons/groups/gursoy_lab/cwalker/projects/pquant/workflow/data/test_fastqs/5k_random_protein_coding_genes.genes_only.fq"
K=25
THRES=0.00001
N_GENES=5000
N_BATCH=250
MEM="200G"

# Create folders
BASE_FOLDER=$"../out/${DATA}"
KMER_FOLDER="${BASE_FOLDER}/${SLURM_JOB_ID}/kmer"
BFV_FOLDER="${BASE_FOLDER}/${SLURM_JOB_ID}/bfv"
SLURM_OUT_FOLDER="${BASE_FOLDER}/${SLURM_JOB_ID}/slurm_out"
STEP4_FOLDER="${BASE_FOLDER}/${SLURM_JOB_ID}/slurm_out/step4"
# create folders
mkdir -p "${BASE_FOLDER}"
mkdir -p "${BASE_FOLDER}/${SLURM_JOB_ID}"
mkdir -p ${KMER_FOLDER}
mkdir -p ${BFV_FOLDER}
mkdir -p ${SLURM_OUT_FOLDER}
mkdir -p ${STEP4_FOLDER}

# GLOBAL_ARGUMENT="-d ${DATA} -k ${K} --thres ${THRES} --kmerFolder ${KMER_FOLDER} --BFVFolder ${BFV_FOLDER} -b"
GLOBAL_ARGUMENT="-k ${K} --thres ${THRES} --kmerFolder ${KMER_FOLDER} --BFVFolder ${BFV_FOLDER} -b -g ${GENE_PATH} -r ${READ_PATH}"

# Step 1
JOBID_STEP1=$(sbatch --output="${SLURM_OUT_FOLDER}/step1.txt" --mem=${MEM} --parsable run_single_step.sh STEP1 "${GLOBAL_ARGUMENT}" )
echo "STEP1: submitted job ${JOBID_STEP1}"

# Step 2 (waits for Step 1 to finish)
JOBID_STEP2=$(sbatch --output="${SLURM_OUT_FOLDER}/step2.txt" --mem=${MEM} --dependency=afterany:${JOBID_STEP1} --parsable run_single_step.sh STEP2 "${GLOBAL_ARGUMENT}")
echo "STEP2: submitted job ${JOBID_STEP2}"

# Step 3 (waits for Step 2 to finish)
JOBID_STEP3=$(sbatch --output="${SLURM_OUT_FOLDER}/step3.txt" --mem=${MEM} --dependency=afterany:${JOBID_STEP2} --parsable run_single_step.sh STEP3 "${GLOBAL_ARGUMENT}")
echo "STEP3: submitted job ${JOBID_STEP3}"

# Step 4 (waits for Step 3 to finish)
# divide N_GENES by NUM_BATCH
# set GENE_START and GENE_END, and add argument by --gs and --ge
N_GENES_PER_BATCH=$(echo "(${N_GENES} - 1) / ${N_BATCH} + 1" | bc)
STEP4_JOB_IDS=()  # Initialize an array to store Step 4 job IDs

for i in $(seq 1 ${N_BATCH})
do
    GENE_START=$(((${i} - 1) * ${N_GENES_PER_BATCH}))
    GENE_END=$(((${i} * ${N_GENES_PER_BATCH}) - 1))
    if [ ${GENE_END} -ge ${N_GENES} ]
    then
        GENE_END=$((${N_GENES} - 1))
    fi
    echo "gene start: ${GENE_START}, gene end: ${GENE_END}"
    JOBID_STEP4=$(sbatch --output="${STEP4_FOLDER}/step4_${i}.txt" --mem=${MEM} --dependency=afterany:${JOBID_STEP3} --array=${i} --parsable run_single_step.sh STEP4 "${GLOBAL_ARGUMENT} --gs ${GENE_START} --ge ${GENE_END}")
    STEP4_JOB_IDS+=(${JOBID_STEP4})  # Store the job ID in the array
    echo "STEP4 ${i}: submitted job ${JOBID_STEP4}"
done

# Step 5 (waits for all Step 4 to finish)
DEPENDENCY_STEP4=$(IFS=":"; echo "${STEP4_JOB_IDS[*]}")
echo "Dependency for STEP5: ${DEPENDENCY_STEP4}"  # Debug statement
echo "gene start: 0, gene end: ${N_GENES}"
JOBID_STEP5=$(sbatch --output="${SLURM_OUT_FOLDER}/step5.txt" --mem=${MEM} --dependency=afterany:${DEPENDENCY_STEP4} --parsable run_single_step.sh STEP5 "${GLOBAL_ARGUMENT} --gs 0 --ge ${N_GENES}")
echo "STEP5: submitted job ${JOBID_STEP5}"

