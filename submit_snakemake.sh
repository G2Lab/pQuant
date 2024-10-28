#!/bin/bash
#SBATCH --job-name=snakemake_workflow       # Job name
#SBATCH --error=slurm_out/snakemake-%j.err  # Error file for Snakemake logs
#SBATCH --output=slurm_out/snakemake-%j.out # Output file for Snakemake logs
#SBATCH --cpus-per-task=8                   # Number of CPUs for Snakemake process
#SBATCH --mem=16G                           # Memory to request for the Snakemake process
#SBATCH --mail-type=ALL
#SBATCH --mail-user=shong@nygenome.org

# Set job_id to the SLURM job ID
job_id=$SLURM_JOB_ID
OUT_DIR=out_test_241021
mkdir -p ${OUT_DIR}/${job_id}/slurm

# Run Snakemake with SLURM integration, using the slurm-jobscript.sh for each rule submission
snakemake -j 256 --config job_id=${job_id} OUT_DIR=${OUT_DIR} \
  --cluster "sbatch --cpus-per-task={resources.cpus} --mem={resources.mem_mb}M --job-name={rule} \
    --time=24:00:00 --output=${OUT_DIR}/${job_id}/slurm/{rule}-%j.out" \
  --latency-wait 60 --printshellcmds --use-conda