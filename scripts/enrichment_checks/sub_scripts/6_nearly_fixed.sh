#!/bin/bash
module load bedtools
# This ends the first group, the completely fixed.
# Next, the nearly fixed. For these we need the vcfs. 
# This script assumes the vcf files are in the vcfs folder here, split into chromosomes.
chr=$1
# Two groups here
# 1. Anc==Ref, AF>90% in vcf
# Pool of positions from the ancestral states
bedtools intersect -a anc_fasta_combined/anc_is_ref_fasta_chr${chr}.bed -b vcfs/sl5_chr${chr}.vcf.bed -sorted -wa -wb | awk '{if($11>0.9) print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$11}'  > nearly_fixed/anc_is_ref_vcf_90_chr${chr}.bed
# 2. Anc!=Ref, AF<10% in vcf
bedtools intersect -a anc_fasta_combined/anc_not_ref_fasta_chr${chr}.bed -b vcfs/sl5_chr${chr}.vcf.bed -sorted -wa -wb | awk '{if($11<0.1) print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$11}' > nearly_fixed/anc_not_ref_vcf_10_chr${chr}.bed
