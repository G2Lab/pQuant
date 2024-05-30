#!/usr/bin/env python
import gzip
import argparse

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


parser = argparse.ArgumentParser()

parser.add_argument(
    "--reads1",
    help="Path to a gzipped fastq file containing R1 reads.",
    default="",
    required=True,
    type=str,
    dest="reads1",
)

parser.add_argument(
    "--reads2",
    help="Path to a gzipped fastq file containing R2 reads.",
    default="",
    required=True,
    type=str,
    dest="reads2",
)


parser.add_argument(
    "--out",
    help="Path to write the merged output to.",
    default="",
    required=True,
    type=str,
    dest="output_path",
)


def smart_open(file_path):
    if file_path.endswith('.gz'):
        return gzip.open(file_path, "rt")
    else:
        return open(file_path, "r")


def main(args):

    file1_path = args.reads1
    file2_path = args.reads2
    output_path = args.output_path


    # Open the gzipped fastq files
    #with gzip.open(file1_path, 'rt') as file1, gzip.open(file2_path, 'rt') as file2, open(output_path, 'w') as outfile:
    with smart_open(file1_path) as file1, smart_open(file2_path) as file2, open(output_path, 'w') as outfile:
        # Initialize fastq iterators
        fastq_iter1 = SeqIO.parse(file1, "fastq")
        fastq_iter2 = SeqIO.parse(file2, "fastq")
        
        for record1, record2 in zip(fastq_iter1, fastq_iter2):
            # Take reverse complement of the second read
            #record2_seq_rc = record2.seq.reverse_complement()
            #record2_qual_rc = record2.letter_annotations["phred_quality"][::-1]
            
            # Concatenate sequences and quality scores
            #concatenated_sequence = record1.seq + record2_seq_rc
            concatenated_sequence = record1.seq + record2.seq
            #concatenated_quality = record1.letter_annotations["phred_quality"] + record2_qual_rc
            concatenated_quality = record1.letter_annotations["phred_quality"] + record2.letter_annotations["phred_quality"]
            
            # Create a new SeqRecord for the concatenated sequence
            concatenated_record = SeqRecord(
                Seq(str(concatenated_sequence)),
                id=record1.id,
                description=record1.description,
                letter_annotations={"phred_quality": concatenated_quality}
            )
            
            # Write concatenated record to the output file
            SeqIO.write(concatenated_record, outfile, "fastq")

    print(f"Concatenated reads are written to {output_path}")


if __name__ == "__main__":
    args = parser.parse_args()
    main(args)
