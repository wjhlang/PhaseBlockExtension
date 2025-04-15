#!/bin/bash
# -----------------------------------------------------------------------------
# Usage: extract_snp_info.sh <input_vcf_or_vcf.gz> <output_txt>
# Example:
#   bash extract_snp_info.sh phased_variants.vcf.gz snp_info.txt
#
# This script:
#   1) Filters for SNPs that are phased heterozygotes from the input VCF.
#   2) Extracts CHROM, POS, REF, ALT, GT, and PS fields.
#   3) Outputs a tab-delimited text file.
#
# It accepts both uncompressed (.vcf) and compressed (.vcf.gz) VCF files.
# -----------------------------------------------------------------------------

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 input_vcf output_txt"
    exit 1
fi

# Assign command-line arguments to variables
INPUT_VCF="$1"
OUTPUT_TXT="$2"

# Execute the bcftools pipeline:
# - Filter for phased heterozygous SNPs.
# - Extract fields: CHROM, POS, REF, ALT, GT, and PS.
bcftools view -v snps -p -i 'GT="het"' "$INPUT_VCF" \
| bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\t[%PS]\n' \
> "$OUTPUT_TXT"
