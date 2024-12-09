import os
from datetime import datetime

configfile: "config/config.yaml"

job_id = config.get("job_id", datetime.now().strftime("%y%m%d-%H%M%S"))  # Use passed job_id or fallback to timestamp

OUT_DIR=f"{config['OUT_DIR']}/{job_id}"
GLOBAL_ARGS = f"-k {config['K']} --thres {config['THRES']} --kmer_folder {OUT_DIR}/kmer --bfv_folder {OUT_DIR}/bfv -g {config['GENE_PATH']} -r {config['READ_PATH']} -b"

rule all:
    input:
        f"{OUT_DIR}/tmp/.run_all_complete"
    output:
        f"{OUT_DIR}/time_summary.csv"
        f"{OUT_DIR}/data_summary/post_statistics.txt"
    resources:
        mem_mb=100,
        cpus=1
    shell:
        """
        ./build/pquant -t tpm {GLOBAL_ARGS} > {OUT_DIR}/data_summary/post_statistics.txt
        python summary/summary.py {OUT_DIR}
        """

rule step1_generate_kmerTable_cloud:
    input:
        f"{OUT_DIR}/tmp/.job_id_initialized"
    output: 
        f"{OUT_DIR}/tmp/step1_done"
    params:
        exe=config['exe']        
    resources:
        mem_mb=int(config['slurm']['step1']['mem'].strip('G')) * 1000, 
        cpus=config['slurm']['step1']['cpus']
    shell:
        """
        {params.exe} -t STEP1 {GLOBAL_ARGS} > {OUT_DIR}/log/step1_generate_kmerTable_cloud.txt
        touch {output}
        """


rule data_statistics:
    input:
        f"{OUT_DIR}/tmp/step1_done"
    output:
        f"{OUT_DIR}/data_summary/data_statistics.txt"
    resources:
        mem_mb=100,
        cpus=1
    shell:
        """
        mkdir -p {OUT_DIR}/data_summary
        ./build/pquant -t analysis {GLOBAL_ARGS} > {OUT_DIR}/data_summary/data_statistics.txt
        """

rule step2_he_keygen_local:
    input:
        rules.step1_generate_kmerTable_cloud.output,
        rules.data_statistics.output
    output:
        f"{OUT_DIR}/tmp/step2_done"
    params:
        exe=config['exe']        
    resources:
        mem_mb=int(config['slurm']['step2']['mem'].strip('G')) * 1000, 
        cpus=config['slurm']['step2']['cpus']
    shell:
        """
        {params.exe} -t STEP2 {GLOBAL_ARGS} > {OUT_DIR}/log/step2_he_keygen_local.txt
        touch {output}
        """

rule step3_encode_and_encrypt_local:
    input:
        rules.step2_he_keygen_local.output
    output:
        f"{OUT_DIR}/tmp/step3_done"
    params:
        exe=config['exe'],
        step3_target=lambda wildcards: "-t STEP3-sim" if config.get('sr') else "-t STEP3",  # Conditional target for STEP3
         sim_args=lambda wildcards: f"--sim_num {config['sim_num']} --sim_len {config['sim_len']}" if config.get('sr') else ""
    resources:
        mem_mb=int(config['slurm']['step3']['mem'].strip('G')) * 1000, 
        cpus=config['slurm']['step3']['cpus']
    shell:
        """
        {params.exe} {params.step3_target} {GLOBAL_ARGS} {params.sim_args} > {OUT_DIR}/log/step3_encode_and_encrypt_local.txt
        touch {output}
        """

rule step4_compute_matching_batch_cloud:
    input:
        rules.step3_encode_and_encrypt_local.output
    output:
        f"{OUT_DIR}/log/step4/step4_{{batch}}.txt"
    params:
        exe=config['exe'],
        N_BATCH=config['N_BATCH']
    resources:
        mem_mb=int(config['slurm']['step4']['mem'].strip('G')) * 1000, 
        cpus=config['slurm']['step4']['cpus']
    shell:
        """
        mkdir -p {OUT_DIR}/log/step4
        echo "Running batch {wildcards.batch}"
        {params.exe} -t STEP4 {GLOBAL_ARGS} --bt {params.N_BATCH} --bn {wildcards.batch} > {OUT_DIR}/log/step4/step4_{wildcards.batch}.txt
        """

rule step4_compute_matching_all_cloud:
    input:
        expand(f"{OUT_DIR}/log/step4/step4_{{batch}}.txt", batch=range(config['N_BATCH']))
    output:
        f"{OUT_DIR}/tmp/step4_compute_matching_cloud.txt"
    resources:
        mem_mb=100,
        cpus=1
    shell:
        """
        touch {output}
        """

rule step5_decrypt_and_return_gene_vector_local:
    input:
        rules.step4_compute_matching_all_cloud.output
    output:
        f"{OUT_DIR}/tmp/.run_all_complete"
    params:
        exe=config['exe']
    resources:
        mem_mb=int(config['slurm']['step5']['mem'].strip('G')) * 1000, 
        cpus=config['slurm']['step5']['cpus']
    shell:
        """
        mkdir -p {OUT_DIR}/bfv/ctxtout
        {params.exe} -t STEP5 {GLOBAL_ARGS} > {OUT_DIR}/log/step5_decrypt_and_return_gene_vector_local.txt
        touch {output}
        """

rule init:
    output:
        f"{OUT_DIR}/tmp/.job_id_initialized",
        f"{OUT_DIR}/parameters.txt"
    params:
        K=config['K'],
        THRES=config['THRES'],
        GENE_PATH=config['GENE_PATH'],
        READ_PATH=config['READ_PATH'],
        N_BATCH=config['N_BATCH'],
        sim_num =lambda wildcards: f"sim_num: {config['sim_num']}" if config.get('sr') else "",
        sim_len =lambda wildcards: f"sim_len: {config['sim_len']}" if config.get('sr') else ""
    resources:
        mem_mb=100,  # Memory in GB
        cpus=1     # Number of CPUs
    shell:
        """
        mkdir -p {OUT_DIR} && mkdir -p {OUT_DIR}/tmp && mkdir -p {OUT_DIR}/log
        touch {output[0]} && touch {output[1]}
        echo " === Parameters === " > {output[1]}
        echo "job_id: {job_id}" >> {output[1]}
        echo "K: {params.K}" >> {output[1]}
        echo "THRES: {params.THRES}" >> {output[1]}
        echo "GENE_PATH: {params.GENE_PATH}" >> {output[1]}
        echo "READ_PATH: {params.READ_PATH}" >> {output[1]}
        echo "N_BATCH: {params.N_BATCH}" >> {output[1]}
        echo "{params.sim_num}" >> {output[1]}
        echo "{params.sim_len}" >> {output[1]}
        echo " ================== " >> {output[1]}
        """