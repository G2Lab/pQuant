#!/bin/bash
#SBATCH --job-name=pq_s1    # Job name
#SBATCH --partition=pe2   # Specify the partition
#SBATCH --nodes=1                 # Number of nodes
#SBATCH --ntasks-per-node=1       # Number of tasks (MPI processes) per node
#SBATCH --cpus-per-task=1         # Number of CPU cores per task
#SBATCH --mail-user=shong@nygenome.org   # Specify your email address

# Navigate to the directory containing your executable
cd /gpfs/commons/groups/gursoy_lab/shong/Github/pQuant_enc/build

./pquant -t $1 $2