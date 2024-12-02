#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 10x_PS_bed PS_0907_bed merged_filtered_bed"
    exit 1
fi

# Assign command-line arguments to variables
PS_10X="$1"
PS_0907="$2"
MERGED_FILTERED_BED="$3"

# Process the pipeline
cat \
    <(awk '{print $0 "\t1\tcontig"}' "$PS_0907") \
    <(bedtools subtract -a "$PS_10X" -b "$PS_0907" -f 1.0 | awk '{print $0 "\t1\t10x"}') \
| sort -k1,1 -k2,2n \
| bedtools merge -c 4,5 -o sum,collapse -i - \
| awk '$4 > 1 && $5 !~ /^(10x,)*10x$/' \
> "$MERGED_FILTERED_BED"
