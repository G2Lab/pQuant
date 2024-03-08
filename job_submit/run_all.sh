#!/bin/bash
#SBATCH --job-name=pq_all    # Job name
#SBATCH --partition=pe2   # Specify the partition
#SBATCH --nodes=1                 # Number of nodes
#SBATCH --ntasks-per-node=1       # Number of tasks (MPI processes) per node
#SBATCH --cpus-per-task=1         # Number of CPU cores per task
#SBATCH --mem=40G                  # Memory per node (adjust as needed)
#SBATCH --mail-user=shong@nygenome.org   # Specify your email address

# Navigate to the directory containing your executable
cd $(pwd)

echo "SLURM_JOB_ID: ${SLURM_JOB_ID}"
# Define paths and arguments
# Initialize variables with default values
K=15
THRES=0.001
N_GENES=5000
N_BATCH=250
MEM="200G"
GENE_PATH="/gpfs/commons/groups/gursoy_lab/cwalker/projects/pquant/workflow/data/reference/pquant/5k_random_protein_coding_genes.combined_exons.exons.fa"
READ_PATH="/gpfs/commons/groups/gursoy_lab/cwalker/projects/pquant/workflow/data/test_fastqs/5k_random_protein_coding_genes.genes_only.fq"
DATA=""
OUT_DIR="../out/021624"

# Parse command line arguments
while getopts ":k:g:r:t:n:b:m:d:o:" opt; do
  case ${opt} in
    k )
      K=$OPTARG
      ;;
    g )
      GENE_PATH=$OPTARG
      ;;
    r )
      READ_PATH=$OPTARG
      ;;
    t )
      THRES=$OPTARG
      ;;
    n )
      N_GENES=$OPTARG
      ;;
    b )
      N_BATCH=$OPTARG
      ;;
    m )
      MEM=$OPTARG
      ;;
    d )
      DATA=$OPTARG
      ;;
    o )
      OUT_DIR=$OPTARG
      ;;    
    \? )
      echo "Invalid option: $OPTARG" 1>&2
      exit 1
      ;;
    : )
      echo "Invalid option: $OPTARG requires an argument" 1>&2
      exit 1
      ;;
  esac
done
shift $((OPTIND -1))

DATA_PATHS=""
# Assign gene and read paths based on the DATA value
if [[ $DATA == "" ]]; then
    DATA_PATHS="-g ${GENE_PATH} -r ${READ_PATH}"
else
    DATA_PATHS="-d ${DATA}"
fi

# Check if either DATA or GENE_PATH and READ_PATH are provided
if ([ -z $DATA ] && ([ -z $GENE_PATH ] || [ -z $READ_PATH ])); then
    echo "Usage: $0 [-d <DATA> | (-g <GENE_PATH> -r <READ_PATH>)] -k <K> -t <THRES> -n <N_GENES> -b <N_BATCH> -m <MEM> -out <OUT_DIR>"
    echo "either DATA or GENE_PATH and READ_PATH must be provided"
    exit 1
fi

# Print out the parameters
echo "DATA_PATHS : $DATA_PATHS"
echo "K: $K"
echo "THRES: $THRES"
echo "N_GENES: $N_GENES"
echo "N_BATCH: $N_BATCH"
echo "MEM: $MEM"
echo "OUT_DIR: $OUT_DIR"


# Create folders
BASE_FOLDER="${OUT_DIR}"
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

GLOBAL_ARGUMENT="-k ${K} --thres ${THRES} --kmer_folder ${KMER_FOLDER} --bfv_folder ${BFV_FOLDER} -b -j ${DATA_PATHS}"

echo "GLOBAL_ARGUMENT: ${GLOBAL_ARGUMENT}"

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

