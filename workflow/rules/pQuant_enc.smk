# Functions
def batch_custom_range_string(start, stop, step):
    values = list(range(start, stop, step))
    if values[-1] != stop:
        values.append(stop)#
    return_strings = []
    for i in values:
        return_strings.append(str(i-step)+"_"+str(i-1))
    return return_strings
    

# Global variables
total_genes = count_unique_ids()
#total_genes = 400
batch_size = config["pquant_batch_size"]


# Rules
rule pQuant_enc_step_one:
    input:
        pquant_fasta="data/reference/pquant/{reference_prefix}.exons.fa"
    output:
        pquant_kmer_table="data/reference/pQuant_enc/{reference_prefix}-k_{k}-entropy_{entropy}/kmer/kmertable_{reference_prefix}-k_{k}-entropy_{entropy}.bin"
    log:
        out="logs/pQuant_enc_step_1/{reference_prefix}.k_{k}_entropy_{entropy}.out",
        err="logs/pQuant_enc_step_1/{reference_prefix}.k_{k}_entropy_{entropy}.err",
    threads:
        8
    resources:
        cpus=8,
        mem_mb=lambda _, attempt: 200000 + ((attempt - 1) * 50000),
        partition="pe2"
    conda:
        "pq"
    shell:
        """
        exec 1> {log.out};
        exec 2> {log.err};

        module load cmake/3.16.3
        module load gcc/9.2.0
        module load clang

        tools/pQuant_enc/build/pquant \
        -k {wildcards.k} \
        -e {wildcards.entropy} \
        -g {input.pquant_fasta} \
        --kmer_folder data/reference/pQuant_enc/{config[reference_prefix]}-k_{wildcards.k}-entropy_{wildcards.entropy}/kmer \
        --bfv_folder data/reference/pQuant_enc/{config[reference_prefix]}-k_{wildcards.k}-entropy_{wildcards.entropy}/bfv \
        -d {config[reference_prefix]}-k_{wildcards.k}-entropy_{wildcards.entropy} \
        -t STEP1
        """


rule pQuant_enc_step_two:
    input:
        pquant_fasta="data/reference/pquant/{reference_prefix}.exons.fa",
        pquant_kmer_table="data/reference/pQuant_enc/{reference_prefix}-k_{k}-entropy_{entropy}/kmer/kmertable_{reference_prefix}-k_{k}-entropy_{entropy}.bin"
    output:
        public_key="data/reference/pQuant_enc/{reference_prefix}-k_{k}-entropy_{entropy}/bfv/key-public.txt",
        private_key="data/reference/pQuant_enc/{reference_prefix}-k_{k}-entropy_{entropy}/bfv/key-private.txt",
    log:
        out="logs/pQuant_enc_step_2/{reference_prefix}.k_{k}_entropy_{entropy}.out",
        err="logs/pQuant_enc_step_2/{reference_prefix}.k_{k}_entropy_{entropy}.err",
    threads:
        8
    resources:
        cpus=8,
        mem_mb=lambda _, attempt: 100000 + ((attempt - 1) * 50000),
        partition="pe2"
    conda:
        "pq"
    shell:
        """
        exec 1> {log.out};
        exec 2> {log.err};

        module load cmake/3.16.3
        module load gcc/9.2.0
        module load clang

        tools/pQuant_enc/build/pquant \
        -k {wildcards.k} \
        -e {wildcards.entropy} \
        -r " " \
        -g {input.pquant_fasta} \
        --kmer_folder data/reference/pQuant_enc/{config[reference_prefix]}-k_{wildcards.k}-entropy_{wildcards.entropy}/kmer \
        --bfv_folder data/reference/pQuant_enc/{config[reference_prefix]}-k_{wildcards.k}-entropy_{wildcards.entropy}/bfv \
        -d {config[reference_prefix]}-k_{wildcards.k}-entropy_{wildcards.entropy} \
        -t STEP2
        """


rule pQuant_enc_step_three:
    input:
        reads=get_pquant_paired_concat_fq,
        pquant_fasta="data/reference/pquant/{reference_prefix}.exons.fa",
        public_key="data/reference/pQuant_enc/{reference_prefix}-k_{k}-entropy_{entropy}/bfv/key-public.txt",
        private_key="data/reference/pQuant_enc/{reference_prefix}-k_{k}-entropy_{entropy}/bfv/key-private.txt",
    output:
        ctxt_read="data/reference/pQuant_enc/{reference_prefix}-k_{k}-entropy_{entropy}/bfv/{sample}/ctxt_read/ct_0.txt",
    log:
        out="logs/pQuant_enc_step_3/{reference_prefix}_{sample}_k_{k}_entropy_{entropy}.out",
        err="logs/pQuant_enc_step_3/{reference_prefix}_{sample}_k_{k}_entropy_{entropy}.err",
    threads:
        8
    resources:
        cpus=8,
        mem_mb=lambda _, attempt: 200000 + ((attempt - 1) * 50000),
        partition="pe2"
    conda:
        "pq"
    shell:
        """
        exec 1> {log.out};
        exec 2> {log.err};

        module load cmake/3.16.3
        module load gcc/9.2.0
        module load clang

        tools/pQuant_enc/build/pquant \
        -k {wildcards.k} \
        -e {wildcards.entropy} \
        -r {input.reads} \
        -g {input.pquant_fasta} \
        --kmer_folder data/reference/pQuant_enc/{config[reference_prefix]}-k_{wildcards.k}-entropy_{wildcards.entropy}/kmer \
        --bfv_folder data/reference/pQuant_enc/{config[reference_prefix]}-k_{wildcards.k}-entropy_{wildcards.entropy}/bfv \
        --ctxt_read_folder data/reference/pQuant_enc/{config[reference_prefix]}-k_{wildcards.k}-entropy_{wildcards.entropy}/bfv/{wildcards.sample}/ctxt_read \
        --ctxt_out_folder data/reference/pQuant_enc/{config[reference_prefix]}-k_{wildcards.k}-entropy_{wildcards.entropy}/bfv/{wildcards.sample}/ctxt_out \
        -d {config[reference_prefix]}-k_{wildcards.k}-entropy_{wildcards.entropy} \
        -t STEP3
        """


adjusted_batches = [int(batch) - int(config["pquant_batch_size"]) for batch in custom_range(batch_size, total_genes, batch_size)]


rule create_batch_files_for_pQuant_enc_step_four:
    input:
        ctxt_read="data/reference/pQuant_enc/{reference_prefix}-k_{k}-entropy_{entropy}/bfv/{sample}/ctxt_read/ct_0.txt" 
    output:
        expand(
            "data/reference/pQuant_enc/{{reference_prefix}}-k_{{k}}-entropy_{{entropy}}/bfv/{{sample}}/batches/batch_{batch}.txt",
            batch=batch_custom_range_string(batch_size, total_genes, batch_size)
        )
    shell:
        """
        touch {output} 
        """


def get_read_counts(wildcards):
    print("wildcards.batch:", wildcards.batch)
    try:
        batch = int(wildcards.batch)
    except ValueError:
        raise ValueError("Invalid value for batch wildcard. It should be an integer.")
    return [batch - i for i in range(5)]


rule pQuant_enc_batched_step_four:
    input:
        pquant_fasta="data/reference/pquant/{reference_prefix}.exons.fa",
        batch_file="data/reference/pQuant_enc/{reference_prefix}-k_{k}-entropy_{entropy}/bfv/{sample}/batches/batch_{gene_start}_{gene_end}.txt",
        # from step 1
        pquant_kmer_table="data/reference/pQuant_enc/{reference_prefix}-k_{k}-entropy_{entropy}/kmer/kmertable_{reference_prefix}-k_{k}-entropy_{entropy}.bin",
        # from step 2
        public_key="data/reference/pQuant_enc/{reference_prefix}-k_{k}-entropy_{entropy}/bfv/key-public.txt",
        private_key="data/reference/pQuant_enc/{reference_prefix}-k_{k}-entropy_{entropy}/bfv/key-private.txt",
        # from step 3
        ctxt_read="data/reference/pQuant_enc/{reference_prefix}-k_{k}-entropy_{entropy}/bfv/{sample}/ctxt_read/ct_0.txt",
    output:
        "data/reference/pQuant_enc/{reference_prefix}-k_{k}-entropy_{entropy}/bfv/{sample}/batch_status/batch_{gene_start}_{gene_end}.done"
    log:
        out="logs/pQuant_enc_step_four/reference_{reference_prefix}.sample_{sample}.k_{k}.entropy_{entropy}.batch_{gene_start}_{gene_end}.out",
        err="logs/pQuant_enc_step_four/reference_{reference_prefix}.sample_{sample}.k_{k}.entropy_{entropy}.batch_{gene_start}_{gene_end}.err"
    threads:
        8
    resources:
        cpus=8,
        mem_mb=lambda _, attempt: 105000 + ((attempt - 1) * 10000),
        partition="pe2"
    conda:
        "pq"
    shell:
        """
        exec 1> {log.out};
        exec 2> {log.err};

        module load cmake/3.16.3
        module load gcc/9.2.0
        module load clang

        mkdir -p data/reference/pQuant_enc/{config[reference_prefix]}-k_{wildcards.k}-entropy_{wildcards.entropy}/bfv/{wildcards.sample}/ctxt_out

        tools/pQuant_enc/build/pquant \
        -k {wildcards.k} \
        -e {wildcards.entropy} \
        -r " " \
        -g {input.pquant_fasta} \
        --kmer_folder data/reference/pQuant_enc/{config[reference_prefix]}-k_{wildcards.k}-entropy_{wildcards.entropy}/kmer \
        --bfv_folder data/reference/pQuant_enc/{config[reference_prefix]}-k_{wildcards.k}-entropy_{wildcards.entropy}/bfv \
        --ctxt_read_folder data/reference/pQuant_enc/{config[reference_prefix]}-k_{wildcards.k}-entropy_{wildcards.entropy}/bfv/{wildcards.sample}/ctxt_read \
        --ctxt_out_folder data/reference/pQuant_enc/{config[reference_prefix]}-k_{wildcards.k}-entropy_{wildcards.entropy}/bfv/{wildcards.sample}/ctxt_out \
        -d {config[reference_prefix]}-k_{wildcards.k}-entropy_{wildcards.entropy} \
        --gs {wildcards.gene_start} \
        --ge {wildcards.gene_end} \
        -t STEP4

        touch {output}
        """


rule pQuant_enc_batched_step_five:
    input:
        pquant_fasta="data/reference/pquant/{reference_prefix}.exons.fa",
        # from step 1
        pquant_kmer_table="data/reference/pQuant_enc/{reference_prefix}-k_{k}-entropy_{entropy}/kmer/kmertable_{reference_prefix}-k_{k}-entropy_{entropy}.bin",
        # from step 2
        public_key="data/reference/pQuant_enc/{reference_prefix}-k_{k}-entropy_{entropy}/bfv/key-public.txt",
        private_key="data/reference/pQuant_enc/{reference_prefix}-k_{k}-entropy_{entropy}/bfv/key-private.txt",
        # from step 3
        ctxt_read="data/reference/pQuant_enc/{reference_prefix}-k_{k}-entropy_{entropy}/bfv/{sample}/ctxt_read/ct_0.txt",
        # from step 4
        batch_statuses=expand(
            "data/reference/pQuant_enc/{{reference_prefix}}-k_{{k}}-entropy_{{entropy}}/bfv/{{sample}}/batch_status/batch_{batch}.done",
            batch=batch_custom_range_string(batch_size, total_genes, batch_size)
        )
    output:
        "data/reference/pQuant_enc/{reference_prefix}-k_{k}-entropy_{entropy}/bfv/{sample}/gene_counts.info"
    log:
        out="logs/pQuant_enc_step_5/reference_{reference_prefix}.sample_{sample}.k_{k}.entropy_{entropy}.out",
        err="logs/pQuant_enc_step_5/reference_{reference_prefix}.sample_{sample}.k_{k}.entropy_{entropy}.err"
    threads:
        8
    resources:
        cpus=8,
        mem_mb=lambda _, attempt: 150000 + ((attempt - 1) * 50000),
        partition="pe2"
    conda:
        "pq"
    shell:
        """
        exec 1> {log.out};
        exec 2> {log.err};

        module load cmake/3.16.3
        module load gcc/9.2.0
        module load clang

        tools/pQuant_enc/build/pquant \
        -k {wildcards.k} \
        -e {wildcards.entropy} \
        -r " " \
        -g {input.pquant_fasta} \
        --kmer_folder data/reference/pQuant_enc/{config[reference_prefix]}-k_{wildcards.k}-entropy_{wildcards.entropy}/kmer \
        --bfv_folder data/reference/pQuant_enc/{config[reference_prefix]}-k_{wildcards.k}-entropy_{wildcards.entropy}/bfv \
        --ctxt_read_folder data/reference/pQuant_enc/{config[reference_prefix]}-k_{wildcards.k}-entropy_{wildcards.entropy}/bfv/{wildcards.sample}/ctxt_read \
        --ctxt_out_folder data/reference/pQuant_enc/{config[reference_prefix]}-k_{wildcards.k}-entropy_{wildcards.entropy}/bfv/{wildcards.sample}/ctxt_out \
        -d {config[reference_prefix]}-k_{wildcards.k}-entropy_{wildcards.entropy} \
        -t STEP5 > {output}
        """

