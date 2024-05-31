import pandas as pd
from glob import glob


##############################################
### Functions for parsing the config files ###
##############################################
samples = (
    pd.read_csv(config["samples"], sep=",", dtype={"sample_name": str, "replicate": str})
    .set_index(["sample_name", "replicate"], drop=False)
    .sort_index()
)

samples_paired = (
    pd.read_csv(config["samples_paired"], sep=",", dtype={"sample_name": str, "replicate": str})
    .set_index(["sample_name", "replicate"], drop=False)
    .sort_index()
)

# code for working with initial pquant testing data, e.g. 5k randomly sampled # protein coding genes
test_samples = (
    pd.read_csv(config["test_samples"], sep=",", dtype={"sample_name": str, "replicate": str})
    .set_index(["sample_name", "replicate"], drop=False)
    .sort_index()
)

# code for working with initial paired end pquant testing data, e.g. 5k randomly sampled # protein coding genes
test_samples_paired = (
    pd.read_csv(config["test_samples_paired"], sep=",", dtype={"sample_name": str, "replicate": str})
    .set_index(["sample_name", "replicate"], drop=False)
    .sort_index()
)

benchmark_samples = (
    pd.read_csv(config["benchmark_samples"], sep=",", dtype={"sample_name": str, "n_reads": str})
    .set_index(["sample_name", "n_reads"], drop=False)
    .sort_index()
)


#############################################
### Functions for parsing the fastq files ###
#############################################
def star_get_subsample_fq(wildcards):
    samp = samples.loc[(wildcards.sample, wildcards.replicate)]
    if samp.paired_end.iloc[0] == 0:
        return samp.fq1
    else:
        return f"{samp.fq1} {samp.fq2}"


def test_star_get_subsample_fq(wildcards):
    samp = test_samples.loc[(wildcards.sample)]
    if int(samp.paired_end.iloc[0]) == 0:
        return samp.fq1
    else:
        return f"{samp.fq1} {samp.fq2}"


# functions to fetch test samples with paired end data
def test_star_get_subsample_fq_paired(wildcards):
    samp = test_samples_paired.loc[(wildcards.sample)]
    print(samp.fq1.item())
    print(samp.fq2.item())
    return f"{samp.fq1} {samp.fq2}"


def test_star_get_subsample_fq_r1(wildcards):
    samp = test_samples_paired.loc[(wildcards.sample)]
    return samp.fq1


def test_star_get_subsample_fq_r2(wildcards):
    samp = test_samples_paired.loc[(wildcards.sample)]
    return samp.fq2


def test_get_pquant_paired_concat_fq(wildcards):
    samp = test_samples_paired.loc[(wildcards.sample)]
    return samp.fqConcat


# functions to fetch samples with paired end data
def get_sample_fq_r1(wildcards):
    samp = samples_paired.loc[(wildcards.sample)]
    return samp.fq1


def get_sample_fq_r2(wildcards):
    samp = samples_paired.loc[(wildcards.sample)]
    return samp.fq2


def get_pquant_paired_concat_fq(wildcards):
    samp = samples_paired.loc[(wildcards.sample)]
    #print(samp)
    return samp.fqConcat


def get_sample_fq(wildcards):
    samp = samples.loc[(wildcards.sample)]
    if int(samp.paired_end.iloc[0]) == 0:
        return samp.fq1
    else:
        return f"{samp.fq1} {samp.fq2}"


def get_sample_name(wildcards):
    return wildcards.sample


def get_sample_benchmark_fq(wildcards):
    samp = benchmark_samples.loc[(wildcards.sample, wildcards.subset_reads)]
    if int(samp.paired_end) == 0:
        return samp.fq1
    else:
        return f"{samp.fq1} {samp.fq2}"

def get_benchmark_paired_concat_fq(wildcards):
    samp = benchmark_samples.loc[(wildcards.sample, wildcards.subset_reads)]
    return samp.fq1


def get_benchmark_subset_gene_list(wildcards):
    gene_list_file = "data/reference/gencode/" + config["reference_prefix"] + "_" + wildcards.gene_subset + ".gene_names.txt"
    return gene_list_file


### pQuant enc ###

# Step 1 input

# Step 4 input
def count_unique_ids():
    ref_prefix = config["reference_prefix"]
    ref_fasta=f"data/reference/pquant/{ref_prefix}.exons.fa"
    unique_ids = set()
    with open(ref_fasta, "r") as f:
        for line in f:
            if line.startswith(">"):
                id = line.strip().split(":")[0][1:]
                unique_ids.add(id)
    return len(unique_ids)



def get_unique_gene_indices_from_ref():
    ref_prefix = config["reference_prefix"]
    ref_fasta=f"data/reference/pquant/{ref_prefix}.exons.fa"
    unique_ids = set()
    with open(ref_fasta, "r") as f:
        for line in f:
            if line.startswith(">"):
                id = line.strip().split(":")[0][1:]
                unique_ids.add(id)
    return [i for i in range(len(unique_ids))]


def generate_gene_index_dictionary(items, sublist_size):
    dictionary = {}
    for i in range(0, len(items), sublist_size):
        sublist = items[i:i+sublist_size]
        dictionary[i//sublist_size] = sublist
    return dictionary


sample_indices = get_unique_gene_indices_from_ref()
pquant_groups = generate_gene_index_dictionary(
    sample_indices,
    50
)
inv_pquant_groups = {}
for x in sample_indices:
    for k in pquant_groups:
        if x in pquant_groups[k]:
            inv_pquant_groups[x] = k


def pQuant_enc_step_four_batch_ids():
    gene_count = count_unique_ids 
    ids = list(range(0, gene_count))
    for i in range(0, len(ids), 50):
        first_id = i
        last_id = min(i + 49, count - 1)
        yield first_id, last_id
        #yield f"batch_{i // 50}.", ids[i:i+50]


def batch_ids():
    count = count_unique_ids()
    for i in range(0, count, 50):
        first_id = i
        last_id = min(i + 49, count - 1)  # Ensure the last batch doesn't exceed the total count
        yield first_id, last_id


#######################################################
### Functions for defining the desired output files ###
#######################################################
def test_final_output():
    targets = []
    for samp in list(set(test_samples_paired["sample_name"])):
        targets.append(f"test_results/STAR_align/{samp}/ReadsPerGene.out.tab")
        targets.append(f"test_results/kallisto/{samp}/abundance.tsv")
        for kmer in [20,25,30,35,40]:
            for entropy in [-1, 0.1, 0.01, 0.001, 0.0001]:
                targets.append(f"test_results/pquant_kmer/{samp}/{samp}.k_{kmer}.entropy_{entropy}.counts.tsv")
    return targets


def test_paired_final_output():
    targets = []
    for samp in list(set(test_samples_paired["sample_name"])):
        print("sample:", samp)
        targets.append(f"test_results/STAR_align_paired/{samp}/ReadsPerGene.out.tab")
        targets.append(f"test_results/kallisto_paired/{samp}/abundance.tsv")
        for kmer in [20,30]:
            for entropy in [0.00001]:
                targets.append(f"test_results/pquant_kmer_sum_paired_r1_only_WITH_RC/{samp}/{samp}.k_{kmer}_entropy_{entropy}.counts.tsv")
                targets.append(f"test_results/pquant_kmer_sum_paired_WITH_RC/{samp}/{samp}.k_{kmer}_entropy_{entropy}.counts.tsv")
                targets.append(f"test_results/pquant_kmer_sum_paired_r1_only/{samp}/{samp}.k_{kmer}_entropy_{entropy}.counts.tsv")
                targets.append(f"test_results/pquant_kmer_sum_paired/{samp}/{samp}.k_{kmer}_entropy_{entropy}.counts.tsv")
    return targets


def benchmark_output():
    targets = []
    for samp in list(set(benchmark_samples["sample_name"])):
        #for subset_reads in range(200000, 5000001, 200000):
            #targets.append(f"benchmark_results/kallisto/{samp}/{samp}_n_{subset_reads}/abundance.tsv")
            #targets.append(f"benchmark_results/STAR_align/{samp}/{samp}_n_{subset_reads}/ReadsPerGene.out.tab")
        for kmer in [20]:
            #for entropy in [-1, 0.1, 0.001]:
            for entropy in [1e-09]:
                for subset_reads in range(5000000, 50000001, 5000000):
                #for subset_reads in [5000000]:
                    targets.append(f"benchmark_results/pquant_kmer/{samp}/{samp}.k_{kmer}_entropy_{entropy}.n_{subset_reads}.counts.tsv")
    return targets


def final_output():
    targets = []
    for samp in list(set(samples["sample_name"])):
        targets.append(f"results/STAR_align/{samp}/ReadsPerGene.out.tab")
        targets.append(f"results/kallisto/{samp}/abundance.tsv")
        #for kmer in [20,25,30]:
        for kmer in [20, 30]:
            #for entropy in [-1, 0.1, 0.001]:
            for entropy in [0.00001]:
                #targets.append(f"results/pquant_kmer/{samp}/{samp}_k_{kmer}_entropy_{entropy}.counts.tsv")
                #targets.append(f"results/pquant_kmer_binary/{samp}/{samp}_k_{kmer}_entropy_{entropy}.counts.tsv")
                for sumthresh in [2,3,4]:
                    targets.append(f"results/pquant_kmer_sum_thresh/{samp}/{samp}_k_{kmer}_entropy_{entropy}.sumthresh_{sumthresh}.counts.tsv")
    return targets


def paired_final_output():
    targets = []
    for samp in list(set(samples_paired["sample_name"])):
        targets.append(f"results/STAR_align_paired/{samp}/ReadsPerGene.out.tab")
        targets.append(f"results/kallisto_paired/{samp}/abundance.tsv")
        for kmer in [15, 16, 17, 18, 19, 20, 25, 30]:
            for entropy in [0.00001]:
                #targets.append(f"results/pquant_kmer_sum_paired_r1_only_WITH_RC/{samp}/{samp}.k_{kmer}.entropy_{entropy}.counts.tsv")
                targets.append(f"results/pquant_kmer_sum_paired_WITH_RC/{samp}/{samp}.k_{kmer}.entropy_{entropy}.counts.tsv")
    return targets


def index_size_benchmark_output():
    target_ref=config["reference_prefix"]
    targets = []
    for gene_subset in range(5000, 60001, 5000):
    #for gene_subset in [5000]:
        # Generate standard kallisto/STAR references, and fasta for pQuant
        targets.append(f"data/reference/kallisto/{target_ref}_{gene_subset}.transcripts.fa")
        targets.append(f"data/reference/star/{target_ref}_{gene_subset}/SA")
        targets.append(f"data/reference/pquant/{target_ref}_{gene_subset}.exons.fa")
        for k in [20]:
            for entropy in [1e-09]:
                # Generate k/entropy-specific encoded reference kmer table
                targets.append(f"data/reference/pquant_encoded/{target_ref}-subset_{gene_subset}-k_{k}-entropy_{entropy}/kmer/kmerTable_{target_ref}-subset_{gene_subset}-k_{k}-entropy_{entropy}.txt")
        
        # Map using kallisto/STAR to get BAM file size

        # Generate encoded/encrpyted pQuant reads
        #targets.append(f"data/reference/pquant_encoded/{target_ref}-subset_{gene_subset}-k_{k}-entropy_{entropy}/bfv/ctxt_read/")
    return targets


def read_size_benchmark_output():
    target_ref=config["reference_prefix"]
    targets = []
    # generate mapped reads in BAM format using kallisto
    for subset_reads in range(5000000, 50000001, 5000000):    
        for samp in list(set(benchmark_samples["sample_name"])):
            targets.append(f"benchmark_results/kallisto/{samp}/{samp}_n_{subset_reads}_ref_subset_60000/abundance.tsv")
            targets.append(f"benchmark_results/STAR_align/{samp}/{samp}_n_{subset_reads}_ref_subset_60000/ReadsPerGene.out.tab")
            # generate encoded/encrypted pQuant reads
            for k in [20]:
                for entropy in [1e-09]:
                    targets.append(f"data/reference/pquant_encoded/{target_ref}-subset_60000-k_{k}-entropy_{entropy}/bfv/ctxt_read-sample_{samp}-subset_reads_{subset_reads}/ct_0.txt")
    return targets



### FINAL PAPER OUTPUT ###

## STAR ##

def STAR_align_paired_output():
    targets = []
    for samp in list(set(samples_paired["sample_name"])):
        targets.append(f"paper_results/STAR/{samp}/ReadsPerGene.out.tab")
    return targets

## kallisto ##

def kallisto_align_paired_output():
    targets = []
    for samp in list(set(samples_paired["sample_name"])):
        targets.append(f"paper_results/kallisto/{samp}/abundance.tsv")
    return targets


## pQuant ##

def pQuant_enc_step_one_output():
    targets = []
    target_ref=config["reference_prefix"]
    for k in [20]:
        for entropy in [1e-09]:
            targets.append(f"data/reference/pQuant_enc/{target_ref}-k_{k}-entropy_{entropy}/kmer/kmertable_{target_ref}-k_{k}-entropy_{entropy}.bin")
    return targets


def pQuant_enc_step_two_output():
    targets = []
    target_ref=config["reference_prefix"]
    for k in [20]:
        for entropy in [1e-09]:
            targets.append(f"data/reference/pQuant_enc/{target_ref}-k_{k}-entropy_{entropy}/bfv/key-public.txt")
            targets.append(f"data/reference/pQuant_enc/{target_ref}-k_{k}-entropy_{entropy}/bfv/key-private.txt")
    return targets


def pQuant_enc_step_three_output():
    targets = []
    target_ref=config["reference_prefix"]
    for samp in list(set(samples_paired["sample_name"])):
        for k in [20]:
            for entropy in [1e-09]:
                targets.append(f"data/reference/pQuant_enc/{target_ref}-k_{k}-entropy_{entropy}/bfv/{samp}/ctxt_read/ct_0.txt")
    return targets


def custom_range(start, stop, step):
    values = list(range(start, stop, step))
    if values[-1] != stop:
        values.append(stop)
    return values


def batch_files_output():
    target_ref = config["reference_prefix"]
    batch_size = config["pquant_batch_size"]
    targets = []
    #total_genes = count_unique_ids()
    total_genes = 400
    for samp in list(set(samples_paired["sample_name"])):
        for k in [20]:
            for entropy in [1e-09]:
                for batch in custom_range(batch_size, total_genes, batch_size):
                    out_f = f"data/reference/pQuant_enc/{target_ref}-k_{k}-entropy_{entropy}/bfv/{samp}/batches/batch_{batch-batch_size}_{batch-1}.txt"
                    targets.append(out_f)
    return targets


def pQuant_enc_step_four_output():
    targets = []
    target_ref=config["reference_prefix"]

    total_genes = count_unique_ids()
    #total_genes = 400

    for samp in list(set(samples_paired["sample_name"])):
        for k in [20]:
            for entropy in [1e-09]:
                for batch in custom_range(batch_size, total_genes, batch_size):
                    out_f = f"data/reference/pQuant_enc/{target_ref}-k_{k}-entropy_{entropy}/bfv/{samp}/batch_status/batch_{batch-batch_size}_{batch-1}.done"
                    targets.append(out_f)
    return targets


def pQuant_enc_step_five_output():
    targets = []
    target_ref=config["reference_prefix"]

    total_genes = count_unique_ids()
    #total_genes = 400

    for samp in list(set(samples_paired["sample_name"])):
        for k in [20]:
            for entropy in [1e-09]:
                out_f = f"data/reference/pQuant_enc/{target_ref}-k_{k}-entropy_{entropy}/bfv/{samp}/gene_counts.info"
                targets.append(out_f)
    return targets


