#!/bin/bash
module load bcftools
# This is a parsing script for the vcfs
# Just make sure what you have in the end is this
# chr/start/end/REF/ALT/AF
# Have it sorted
# Have indels removed

# For lycopersicum
chr=$1
# Take relevant columns from the vcf
bcftools query -f '%CHROM\t%POS\t%POS\t%REF\t%ALT\t%AF\n' <(zcat vcfs/chr${chr}.vcf.gz) | \ 
# Convert to bed format
awk '{print $1"\t"($3-1)"\t"$3"\t"$4"\t"$5"\t"$6}' | \
# Remove indels and masked positions
awk '{if(($4 != "-" && $4 != "N") && ($5 != "-" && $5 != "N")) print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' | \ 
# Sort
sort -k1,1n -k2,2n > vcfs/chr_${chr}.vcf.bed