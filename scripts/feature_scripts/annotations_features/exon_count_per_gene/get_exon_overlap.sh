#!/bin/bash
module load bedtools
module load tabix

gff_file="/home/labs/alevy/omerbar/features/expanded_SL5/gff_features/SL5.gff3"
# Parse gff to bed and genes and split to chromosomes
awk '($3=="exon"){print $1"\t"($4-1)"\t"$5}' ${gff_file}  | awk '{print > "exon_chr_"$1".bed"}' 

# go over chromosomes
for i in {0..12}
do 
    sort -k1,1n -k2,2n exon_chr_${i}.bed > exons_chr_${i}_sorted.bed
    bedtools intersect -a <(zcat ~/features/expanded_SL5/feature_file_generation/SL5_chr_${i}_CAP.bed.gz) -b  exons_chr_${i}_sorted.bed  -c -sorted | \
    awk '{print $1"\t"$2"\t"$3"\t"$5}' |\
    bgzip > exon_overlap_count_feature_chr_${i}.bed.gz &
    rm exons_chr_${i}_sorted.bed exon_chr_${i}.bed
done