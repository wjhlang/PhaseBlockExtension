#!/bin/bash

# -----------------------------------------------------------------------------
# Usage: merge_ps.sh <bridge_ps_bed> <target_ps_bed> <merged_ps_bed>
# Example:
#   bash merge_ps.sh 10x_PS.bed contig_PS.bed merged_ps.bed
#
# This script merges and bridges target plaseblock intervals (contigs) with bridging phaseblock intervals (10x).
# The last awk line shows only intervals where:
#   - sum of "count" is >1, meaning they come from both phaseblocks
#   - the "collapsed" field is not purely "bridge,bridge,..." (i.e., must have target too).
# -----------------------------------------------------------------------------


# Check if the correct number of arguments is provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 bridge_ps_bed target_ps_bed merged_filtered_bed"
    exit 1
fi

# Assign command-line arguments to variables
bridge_ps="$1"
target_ps="$2"
merged_ps="$3"

# Process the pipeline
cat \
    <(awk '{print $0 "\t1\ttarget"}' "$target_ps") \
    <(bedtools subtract -a "$bridge_ps" -b "$target_ps" -f 1.0 | awk '{print $0 "\t1\tbridge"}') \
| sort -k1,1 -k2,2n \
| bedtools merge -c 4,5 -o sum,collapse -i - \
| awk '$4 > 1 && $5 !~ /^(bridge,)*bridge$/' \
> "$merged_ps"
