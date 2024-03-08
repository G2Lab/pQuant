#!/bin/bash
#SBATCH --job-name=table_gen    # Job name
#SBATCH --partition=pe2   # Specify the partition
#SBATCH --nodes=1                 # Number of nodes
#SBATCH --ntasks-per-node=1       # Number of tasks (MPI processes) per node
#SBATCH --cpus-per-task=1         # Number of CPU cores per task
#SBATCH --mem=40G                  # Memory per node (adjust as needed)
#SBATCH --output=table_gen.out    # Standard output and error log
#SBATCH --mail-user=shong@nygenome.org   # Specify your email address

# Navigate to the directory containing your executable
cd /gpfs/commons/groups/gursoy_lab/shong/Github/pQuant_enc/build

./pquant -t table -d 5k -k 15