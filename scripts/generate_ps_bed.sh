#!/bin/bash

# -----------------------------------------------------------------------------
# Usage: Generate_PS_bed.sh <input_vcf> <output_bed>
# Example:
#   bash Generate_PS_bed.sh heter.SNP.sorted.phasetag.0907.vcf.gz contig_PS.bed
#
# This script:
#   1) Filters for SNPs that are phased heterozygotes
#   2) Extracts CHROM, POS, and PS
#   3) Groups consecutive positions with the same PS into a single BED interval
# -----------------------------------------------------------------------------


# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 input_vcf_file output_bed_file"
    exit 1
fi

# Assign command-line arguments to variables
INPUT_VCF="$1"
OUTPUT_BED="$2"

# Run the pipeline
bcftools view -v snps -p -i 'GT="het"' "$INPUT_VCF" \
| bcftools query -f '%CHROM\t%POS\t[%PS]\n' \
| awk '$3 != "." {print $0}' \
| sort -k1,1 -k3,3n -k2,2n \
| awk 'BEGIN {OFS="\t"}
{
    if ($1 != prev_chrom || $3 != prev_ps) {
        if (NR > 1) print prev_chrom, ps_start, prev_pos
        prev_chrom = $1
        prev_ps = $3
        ps_start = $3
    }
    prev_pos = $2
}
END {
    if (NR > 0) print prev_chrom, ps_start, prev_pos
}' > "$OUTPUT_BED"
