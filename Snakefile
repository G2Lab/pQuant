import os
from datetime import datetime

configfile: "config/config.yaml"

job_id = config.get("job_id", datetime.now().strftime("%y%m%d-%H%M%S"))  # Use passed job_id or fallback to timestamp

GLOBAL_ARGS = f"-k {config['K']} --thres {config['THRES']} --kmer_folder out/{job_id}/kmer --bfv_folder out/{job_id}/bfv -g {config['GENE_PATH']} -r {config['READ_PATH']} -b"

rule all:
    input:
        f"out/{job_id}/tmp/.run_all_complete"
    output:
        f"out/{job_id}/time_summary.csv"
    resources:
        mem_mb=100,
        cpus=1
    shell:
        """
        python summary/summary.py {job_id}
        """

rule step1_generate_kmerTable_cloud:
    input:
        f"out/{job_id}/tmp/.job_id_initialized"
    output: 
        f"out/{job_id}/tmp/step1_done"
    params:
        exe=config['exe']        
    resources:
        mem_mb=int(config['slurm']['step1']['mem'].strip('G')) * 1000, 
        cpus=config['slurm']['step1']['cpus']
    shell:
        """
        {params.exe} -t STEP1 {GLOBAL_ARGS} > out/{job_id}/step1_generate_kmerTable_cloud.txt
        touch {output}
        """

rule step2_he_keygen_local:
    input:
        rules.step1_generate_kmerTable_cloud.output
    output:
        f"out/{job_id}/tmp/step2_done"
    params:
        exe=config['exe']        
    resources:
        mem_mb=int(config['slurm']['step2']['mem'].strip('G')) * 1000, 
        cpus=config['slurm']['step2']['cpus']
    shell:
        """
        {params.exe} -t STEP2 {GLOBAL_ARGS} > out/{job_id}/step2_he_keygen_local.txt
        touch {output}
        """

rule step3_encode_and_encrypt_local:
    input:
        rules.step2_he_keygen_local.output
    output:
        f"out/{job_id}/tmp/step3_done"
    params:
        exe=config['exe']        
    resources:
        mem_mb=int(config['slurm']['step3']['mem'].strip('G')) * 1000, 
        cpus=config['slurm']['step3']['cpus']
    shell:
        """
        {params.exe} -t STEP3 {GLOBAL_ARGS} > out/{job_id}/step3_encode_and_encrypt_local.txt
        touch {output}
        """

rule step4_compute_matching_batch_cloud:
    input:
        rules.step3_encode_and_encrypt_local.output
    output:
        f"out/{job_id}/step4/step4_{{batch}}.txt"
    params:
        exe=config['exe'],
        N_BATCH=config['N_BATCH']
    resources:
        mem_mb=int(config['slurm']['step4']['mem'].strip('G')) * 1000, 
        cpus=config['slurm']['step4']['cpus']
    shell:
        """
        mkdir -p out/{job_id}/step4
        echo "Running batch {wildcards.batch}"
        {params.exe} -t STEP4 {GLOBAL_ARGS} --bt {params.N_BATCH} --bn {wildcards.batch} > out/{job_id}/step4/step4_{wildcards.batch}.txt
        """

rule step4_compute_matching_all_cloud:
    input:
        expand(f"out/{job_id}/step4/step4_{{batch}}.txt", batch=range(config['N_BATCH']))
    output:
        f"out/{job_id}/tmp/step4_compute_matching_cloud.txt"
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
        f"out/{job_id}/tmp/.run_all_complete"
    params:
        exe=config['exe']
    resources:
        mem_mb=int(config['slurm']['step5']['mem'].strip('G')) * 1000, 
        cpus=config['slurm']['step5']['cpus']
    shell:
        """
        mkdir -p out/{job_id}/bfv/ctxtout
        {params.exe} -t STEP5 {GLOBAL_ARGS} > out/{job_id}/step5_decrypt_and_return_gene_vector_local.txt
        touch {output}
        """

rule init:
    output:
        f"out/{job_id}/tmp/.job_id_initialized",
        f"out/{job_id}/parameters.txt"
    params:
        K=config['K'],
        THRES=config['THRES'],
        GENE_PATH=config['GENE_PATH'],
        READ_PATH=config['READ_PATH'],
        N_BATCH=config['N_BATCH']
    resources:
        mem_gb=1,  # Memory in GB
        cpus=1     # Number of CPUs
    shell:
        """
        mkdir -p out/{job_id} && mkdir -p out/{job_id}/tmp
        touch {output[0]} && touch {output[1]}
        echo " === Parameters === " > {output[1]}
        echo "job_id: {job_id}" >> {output[1]}
        echo "K: {params.K}" >> {output[1]}
        echo "THRES: {params.THRES}" >> {output[1]}
        echo "GENE_PATH: {params.GENE_PATH}" >> {output[1]}
        echo "READ_PATH: {params.READ_PATH}" >> {output[1]}
        echo "N_BATCH: {params.N_BATCH}" >> {output[1]}
        echo " ================== " >> {output[1]}
        """
