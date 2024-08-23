#!/bin/bash
module load bedtools
# This ends the second group, the nearly fixed.
# Next, the rare mutations. 
# This script assumes the vcf files are in the vcfs folder here, split into chromosomes.
chr=$1
# Two groups here
# 1. Anc==Ref, AF<1% in vcf
# Pool of positions from the ancestral states
bedtools intersect -a anc_fasta_combined/anc_is_ref_fasta_chr${chr}.bed -b vcfs/sl5_chr${chr}.vcf.bed -sorted -wa -wb | awk '{if($11<0.01) print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$11}'  > rare_mutations/anc_is_ref_vcf_1_chr${chr}.bed
# 2. Anc!=Ref, AF>99% in vcf
bedtools intersect -a anc_fasta_combined/anc_not_ref_fasta_chr${chr}.bed -b vcfs/sl5_chr${chr}.vcf.bed -sorted -wa -wb | awk '{if($11>0.99) print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$11}' > rare_mutations/anc_not_ref_vcf_99_chr${chr}.bed
