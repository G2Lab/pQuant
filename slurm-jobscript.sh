#!/bin/bash
#SBATCH --job-name={rule}                   # Job name
#SBATCH --output=slurm_out/{rule}-%j.out    # Standard output log file
#SBATCH --error=slurm_out/{rule}-%j.err     # Standard error log file
#SBATCH --cpus-per-task={resources.cpus}    # Number of CPUs (specified per rule)
#SBATCH --mem={resources.mem_mb}M           # Memory (specified per rule)
#SBATCH --time=24:00:00                     # Time limit per job

{exec_job}