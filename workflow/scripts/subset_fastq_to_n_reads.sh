#!/bin/bash

# Check if correct number of arguments is provided
if [[ "$#" -ne 1 ]]; then
    echo "Usage: $0 <input_fastq>"
    exit 1
fi

# Assign input arguments to variables
inreads="$1"

# Extract filename without extension
filename=$(basename -- "$inreads")
name_no_ext="${filename%.fastq}"

# Specify the output directory
outdir="../data/encode/subsample/$name_no_ext"
outdir="${outdir%/}"
# Create the output directory if it doesn't exist
mkdir -p "$outdir"

# Loop to subsample reads in 200,000 read increments up to 10,000,000
for nreads in $(seq 30000000 5000000 50000000); do
    outreads="${outdir}/${name_no_ext}.n_${nreads}.fastq"
    # Subsample using seqtk and save to new file
    seqtk sample -s 42 "$inreads" "$nreads" > "$outreads"
    echo "Subsampled ${nreads} reads to: ${outreads}"
done


echo "Done!"
