#!/usr/bin/env python

import argparse
import re
import multiprocessing
import subprocess
import pysam

from itertools import groupby
from Bio.Seq import Seq
from os import path, makedirs


parser = argparse.ArgumentParser()

parser.add_argument(
    "--gene-list",
    help="Path to a new-line-separated list of gene IDs.",
    default="",
    required=True,
    type=str,
    dest="gene_list",
)

parser.add_argument(
    "--gtf",
    help="Path to a GENCODE GTF.",
    default="",
    required=False,
    type=str,
    dest="gtf",
)

parser.add_argument(
    "--genome-fasta",
    help="Path to a reference human genome in FASTA format.",
    default="",
    required=False,
    type=str,
    dest="genome_fasta",
)

parser.add_argument(
    "--transcript-fasta",
    help="Path to a reference human genome in FASTA format.",
    default="",
    required=False,
    type=str,
    dest="transcript_fasta",
)

parser.add_argument(
    "--out-dir",
    help="Base output directory for generated output index/reference files.",
    default="",
    required=False,
    type=str,
    dest="output_dir",
)

parser.add_argument(
    "--out-prefix",
    help="File prefix for generated output index/reference files.",
    default="",
    required=False,
    type=str,
    dest="output_prefix",
)


def parse_gene_list(gene_list_file):
    with open(gene_list_file, "r") as gene_file:
        gene_ids = [i.strip() for i in gene_file.readlines()]
    return gene_ids


def read_genes_to_set(gene_list):
    return set(gene_list)


def process_chunk(args):
    chunk, gene_set = args
    gene_entries = {}
    for line in chunk:
        if not line.startswith("#"):
            entry = re.split(r'\t+', line.strip())
            gene_id_match = re.search(r'gene_id "([^"]+)"', entry[8])
            if gene_id_match:
                gene_id = gene_id_match.group(1)
                if gene_id in gene_set:
                    gene_entries.setdefault(gene_id, []).append(line)
    return [line for entry_lines in gene_entries.values() for line in entry_lines]


def create_subset_gtfs(input_gtf_path, gene_list, output_gtf_path):
    gene_set = read_genes_to_set(gene_list)
    chunk_size = 10000

    with open(input_gtf_path, "r") as input_file:
        with open(output_gtf_path, "w") as output_file:
            pool = multiprocessing.Pool()
            chunks = []
            chunk = []
            for line in input_file:
                chunk.append(line)
                if len(chunk) >= chunk_size:
                    chunks.append(chunk)
                    chunk = []
            if chunk:
                chunks.append(chunk)

            args_list = [(chunk, gene_set) for chunk in chunks]
            filtered_results = pool.map(process_chunk, args_list)

            for filtered_lines in filtered_results:
                for line in filtered_lines:
                    output_file.write(line)
    return


def create_filtered_gtf_by_feature(input_gtf_path, feature, output_gtf_path):
    with open(input_gtf_path, "r") as input_file:
        with open(output_gtf_path, "w") as output_file:
            for line in input_file:
                if not line.startswith("#"):
                    entry = line.rstrip().split("\t")
                    if entry[2] == feature:
                        output_file.write(line)
    return


def gene_id_from_header(header):
    return header.split("|")[1]


def fasta_iter(fasta_name, gene_ids_to_extract):
    gene_ids_set = set(gene_ids_to_extract)
    with open(fasta_name, "r") as fh:
        faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
        for header_lines in faiter:
            header = header_lines.__next__()[1:].strip()
            seq = "".join(s.strip() for s in faiter.__next__())
            if gene_id_from_header(header) in gene_ids_set:
                yield header, seq


def write_subset_transcript_fasta(transcript_fasta, gene_ids_to_extract, output_fasta_path):
    with open(output_fasta_path, "w") as output_file:
        for header, sequence in fasta_iter(transcript_fasta, gene_ids_to_extract):
            output_file.write(f">{header}\n{sequence}\n")


def merge_overlapping_exons(exon_entries):
    merged_entries = []
    exon_entries.sort(key=lambda x: (x["seqname"], x["start"], x["end"]))
    
    if not exon_entries:
        return []

    current_entry_merged_exon_counts = {}

    current_entry = exon_entries[0]
    for entry in exon_entries[1:]:
        if entry["gene_id"] == current_entry["gene_id"] and entry["start"] <= current_entry["end"]:
            current_entry["end"] = max(current_entry["end"], entry["end"])
        else:
            merged_entries.append(current_entry)
            current_entry = entry

    merged_entries.append(current_entry)
    return merged_entries


def parse_gtf(gtf_file):
    exon_entries = []
    gene_exon_counts = {}

    with open(gtf_file, "r") as file:
        for line in file:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) != 9:
                continue
            seqname, source, feature, start, end, score, strand, frame, attribute = parts
            if feature == "exon":
                attributes_dict = dict(
                    item.strip().split() for item in attribute.strip().split(";") if item.strip()
                )
                gene_id = attributes_dict.get("gene_id", "unknown").strip('"')

                if gene_id not in gene_exon_counts:
                    gene_exon_counts[gene_id] = 1
                else:
                    gene_exon_counts[gene_id] += 1

                exon_number = gene_exon_counts[gene_id]
                region_id = f"{gene_id}:exon{exon_number}_{seqname}:{start}-{end}"
                exon_entries.append(
                    {
                        "seqname": seqname,
                        "start": int(start),
                        "end": int(end),
                        "strand": strand,
                        "region_id": region_id,
                        "gene_id": gene_id,
                    }
                )

    return merge_overlapping_exons(exon_entries)


def extract_regions(fasta_file, gtf_entries):
    try:
        fasta = pysam.FastaFile(fasta_file)
    except FileNotFoundError:
        raise FileNotFoundError("FASTA file not found.")

    region_sequences = {}

    for entry in gtf_entries:
        chrom = entry["seqname"]
        start = entry["start"]
        end = entry["end"]
        strand = entry["strand"]
        gene_id = entry["region_id"].split(":")[0]
        exon_number = entry["region_id"].split(":exon")[1].split("_")[0]
        region_id = f"{gene_id}:exon{exon_number}_{chrom}:{start}-{end}_{strand}"

        try:
            # Use pysam to extract the sequence from the reference FASTA file
            sequence = fasta.fetch(chrom, start, end)
            if strand == "-":
                sequence = str(Seq(sequence).reverse_complement())
            else:
                # initial run with 5k protein coding genes required a utf8 decode
                # second run had sequence already as a string so decoding raised an 
                # error. could be due to different versions of pysam but just
                # catching and explicitly making the sequence a string (redundant)
                try:
                    sequence = sequence.decode("utf-8")
                except AttributeError:
                    sequence = str(sequence)
            region_sequences[region_id] = sequence
        except KeyError:
            print(f"Warning: Region {region_id} not found in the reference FASTA file.")
    
    return region_sequences


def create_star_index(genome_dir, genome_fasta_file, sjdb_gtf_file, args):
    cmd = [
        "STAR",
        "--runMode",
        "genomeGenerate",
        "--genomeDir",
        genome_dir,
        "--genomeFastaFiles",
        genome_fasta_file,
        "--sjdbGTFfile",
        sjdb_gtf_file,
        "--outTmpDir",
        f"_STARtmp_{args.output_prefix}"

    ]
    subprocess.run(cmd, check=True)
    return


def create_kallisto_index(kallisto_index_path, input_fasta_path):
    try:
        subprocess.run(
            ["kallisto", "index", "-i", kallisto_index_path, input_fasta_path],
            check=True,
        )
        print("Kallisto index created successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error creating Kallisto index: {e}")


def create_pquant_ref(subset_gtf, genome_fasta, out_fasta_path):
    gtf_dic = parse_gtf(subset_gtf)
    regions = extract_regions(genome_fasta, gtf_dic)
    out_lines = []
    for k, v in regions.items():
        out_lines.append(f">{k}\n{v}")
    with open(out_fasta_path, "w") as out_f:
        out_f.write("\n".join(out_lines))
    return


def main(args):
    # define output directories
    gtf_out_dir = path.join(args.output_dir, "gtf")
    star_out_dir = path.join(args.output_dir, "star")
    kallisto_out_dir = path.join(args.output_dir, "kallisto")
    pquant_out_dir = path.join(args.output_dir, "pquant")

    # 1) parse list of gene_ids
    print("Parsing genes...")
    genes_to_subset = parse_gene_list(args.gene_list)

    # 2) create output gtfs
    print("Subsetting input GTF...")
    makedirs(gtf_out_dir, exist_ok=True)
    out_subset_gtf = path.join(gtf_out_dir, args.output_prefix + ".gtf")
    out_subset_gene_only_gtf = path.join(
        gtf_out_dir, args.output_prefix + ".genes_only.gtf"
    )
    out_subset_exon_only_gtf = path.join(
        gtf_out_dir, args.output_prefix + ".exons_only.gtf"
    )
    create_subset_gtfs(
        args.gtf,
        genes_to_subset,
        out_subset_gtf 
    )
    create_filtered_gtf_by_feature(out_subset_gtf, "gene", out_subset_gene_only_gtf)
    create_filtered_gtf_by_feature(out_subset_gtf, "exon", out_subset_exon_only_gtf)
    # 3) create star index using main output gtf
    print("Creating STAR index...")
    makedirs(star_out_dir, exist_ok=True)
    star_genome_dir = path.join(star_out_dir, args.output_prefix)
    create_star_index(
        genome_dir=star_genome_dir,
        genome_fasta_file=args.genome_fasta,
        sjdb_gtf_file=out_subset_gtf,
        args=args
    )

    # 4) create kallisto index using subset gencode transcript fasta
    print("Creating Kallisto index...")
    makedirs(kallisto_out_dir, exist_ok=True)
    out_kallisto_fasta = path.join(
        kallisto_out_dir, args.output_prefix + ".transcripts.fa"
    )
    out_kallisto_index = path.join(
        kallisto_out_dir, args.output_prefix + ".transcripts.idx"
    )
    write_subset_transcript_fasta(
        transcript_fasta=args.transcript_fasta,
        gene_ids_to_extract=genes_to_subset,
        output_fasta_path=out_kallisto_fasta
    )
    create_kallisto_index(
        input_fasta_path=out_kallisto_fasta,
        kallisto_index_path=out_kallisto_index
    )

    # 5) create pquant reference fasta using exon gtf
    print("Creating pQuant reference...")
    makedirs(pquant_out_dir, exist_ok=True)
    out_pquant_fasta = path.join(
        pquant_out_dir, args.output_prefix + ".exons.fa"
    )
    create_pquant_ref(
        out_subset_gtf,
        args.genome_fasta,
        out_pquant_fasta
    )
    print("Done!")
    return


if __name__ == "__main__":
    args = parser.parse_args()
    main(args)

