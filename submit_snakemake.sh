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
    --output=out/${job_id}/slurm/{rule}-%j.err" \
  --latency-wait 60 --printshellcmds --use-conda