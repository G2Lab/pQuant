#!/usr/bin/env python

from sys import argv


def main():
    counts = {}
    input_order = []
    with open(argv[1], "r") as f:
        f.readline()
        for line in f:
            line = line.split()
            gene_id = line[0].split("|")[1]
            if gene_id not in counts:
                counts[gene_id] = 0
                input_order.append(gene_id)
            count = float(line[-2])
            counts[gene_id] += count
    for gene in input_order:
        print(f"{gene}\t{counts[gene]}")


if __name__ == "__main__":
    main()

