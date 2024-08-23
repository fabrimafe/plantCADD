#!/bin/bash
# A script to expand the .mod files to bed files and make them 0 based
module load tabix

chr=$1
fasta_file="/home/labs/alevy/omerbar/features/expanded_SL5/feature_file_generation/SL5_chr_${chr}.bed.gz"
sort -k1,1n -k2,2n parsed_repeats_chr_${chr}.mod | awk '{print $1"\t"$2-1"\t"$3"\t"$4"\t"$5}' > parsed_repeats_chr_${chr}_sorted.bed
bedtools intersect -a <(zcat ${fasta_file}) -b  parsed_repeats_chr_${chr}_sorted.bed -wa -wb -sorted | awk '{print $1"\t"$2"\t"$3"\t"$8"\t"$9}' > repeats_chr_${chr}.bed
rm parsed_repeats_chr_${chr}.mod parsed_repeats_chr_${chr}_sorted.bed
