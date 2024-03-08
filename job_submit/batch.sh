#!/bin/bash
#SBATCH --job-name=batch
#SBATCH --output=out_%j.txt
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=NONE
#SBATCH --mail-user=shong@nygenome.org

GENE_PATH="/gpfs/commons/groups/gursoy_lab/cwalker/projects/pquant/workflow/data/reference/pquant/5k_random_protein_coding_genes.combined_exons.exons.fa"
READ_PATH="/gpfs/commons/groups/gursoy_lab/cwalker/projects/pquant/workflow/data/test_fastqs/5k_random_protein_coding_genes.genes_only.fq"

# ### to run all dataset!
for ((K=15; K<=35; K+=5))
do
    for thres in 0.0000001 0.00001 0.001 0.1 100
    do 
        # filename="${output_directory}/pQuant_${K}_${thres}.out"
        sbatch run_all.sh -k ${K} -g ${GENE_PATH} -r ${READ_PATH} -t ${thres} -n 5000 -b 250 -m "60G" -o "../out/030724"
    done
done
