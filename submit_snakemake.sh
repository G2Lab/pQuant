#!/bin/bash
#SBATCH --job-name=snakemake_workflow       # Job name
#SBATCH --error=out/slurm_normal_out_241112/snakemake-%j.err  # Error file for Snakemake logs
#SBATCH --output=out/slurm_normal_out_241112/snakemake-%j.out # Output file for Snakemake logs
#SBATCH --mem=1G                           # Memory to request for the Snakemake process
#SBATCH --mail-type=ALL
#SBATCH --mail-user=shong@nygenome.org
#SBATCH --time=100:00:00                    # Time limit for the job

# Set job_id to the SLURM job ID
job_id="${SLURM_JOB_ID}"
OUT_DIR=out
mkdir -p ${OUT_DIR}/${job_id}/slurm

# Run Snakemake with SLURM integration, using the slurm-jobscript.sh for each rule submission
snakemake -j 64 --config job_id=${job_id} OUT_DIR=${OUT_DIR} \
  --cluster "sbatch --cpus-per-task={resources.cpus} --mem={resources.mem_mb}M --job-name={rule} \
    --time=100:00:00 \
    --output=out/${job_id}/slurm/{rule}-%j.out \
    --error=out/${job_id}/slurm/{rule}-%j.err" \
  --latency-wait 60 --printshellcmds --use-conda
