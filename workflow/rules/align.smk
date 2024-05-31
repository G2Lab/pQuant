rule STAR_align_paired:
    """Perform a STAR alignment of the raw fastq files.

    Input:
        raw FASTQ file directory path containing the bulk RNA-seq data for the sample.
    Output:
        alignment SAM file and gene counts.
    """
    input:
        reads1=get_sample_fq_r1,
        reads2=get_sample_fq_r2,
        index="data/reference/star/" + config["star_kallisto_reference_prefix"]
    output:
        aln="paper_results/STAR/{sample}/Aligned.sortedByCoord.out.bam",
        reads_per_gene="paper_results/STAR/{sample}/ReadsPerGene.out.tab"
    log:
        out="logs/STAR_align_paired/{sample}.out",
        err="logs/STAR_align_paired/{sample}.err",
    threads:
        8
    resources:
        mem_mb=64000, cpus=8
    shell:
        """
        exec 1> {log.out};
        exec 2> {log.err};

        STAR \
        --runMode alignReads \
        --genomeDir {input.index} \
        --runThreadN {threads} \
        --readNameSeparator space \
        --outSAMtype BAM SortedByCoordinate \
        --genomeLoad NoSharedMemory \
        --quantMode GeneCounts \
        --outFileNamePrefix paper_results/STAR/{wildcards.sample}/ \
        --readFilesCommand zcat \
        --readFilesIn {input.reads1} {input.reads2}
        """


rule kallisto_align:
    input:
        reads1=get_sample_fq_r1,
        reads2=get_sample_fq_r2,
        index="data/reference/kallisto/" + config["star_kallisto_reference_prefix"] + ".transcripts.idx"
    output:
        transcript_counts="paper_results/kallisto/{sample}/abundance.tsv",
        gene_counts="paper_results/kallisto/{sample}/gene_counts.tsv"
    log:
        out="logs/kallisto_align_paired/{sample}.out",
        err="logs/kallisto_align_paired/{sample}.err",
    threads:
        8
    resources:
        mem_mb=64000, cpus=8
    shell:
        """
        exec 1> {log.out};
        exec 2> {log.err};

        kallisto quant \
        -i {input.index} \
        -o paper_results/kallisto/{wildcards.sample} \
        -b 1 \
        {input.reads1} {input.reads2}

        python scripts/collapse_kallisto_counts.py \
            {output.transcript_counts} \
            paper_results/kallisto/{wildcards.sample}/gene_counts.tsv \
            > {output.gene_counts}
        """
