#!/bin/bash

# Continuing, now we intersect the fasta nucleotides with the ancestral nucleotides
# We run this via a for loop on the command line

module load bedtools
module load tabix

current_chrom=$1
nucleotide_file_name=$2
anc_file_name=$3

mkdir -p completely_fixed
mkdir -p anc_fasta_intersect
#---------------------------OLD-------------------------------
# Make the nucleotide files all capital letters - done in the python parsing of the fasta
#echo "Capitalizing fasta nucs..."
#cat fasta_nucleotides/${nucleotide_file_name}_chr${current_chrom}.bed |awk '{ $4 = toupper($4); print $1"\t"$2"\t"$3"\t"$4}' | sort -k1,1n -k2,2n > completely_fixed/tmp_nuc_chr_${current_chrom}.bed
#---------------------------OLD-------------------------------

# Intersect with ancestral file
echo "Intersecting nucs and ancestral files..."
bedtools intersect -a <(zcat fasta_nucleotides/${nucleotide_file_name}_chr_${current_chrom}.bed) -b <(zcat anc_states/${anc_file_name}_chr_${current_chrom}.bed) -sorted -wa -wb | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$8}' > anc_fasta_intersect/tmp_anc_fasta_int_chr_${current_chrom}.bed

# Clean for indels
#echo "Cleaning indels..."
cat anc_fasta_intersect/tmp_anc_fasta_int_chr_${current_chrom}.bed | awk '{if(($4 != "-" && $4 != "N") && ($5 != "-" && $5 != "N")) print $1"\t"$2"\t"$3"\t"$4"\t"$5}' > anc_fasta_intersect/anc_fasta_int_chr_${current_chrom}.bed

# Sorting out all anc!=ref
awk '{if($4==$5) print $1"\t"$2"\t"$3"\t"$4"\t"$5}' anc_fasta_intersect/anc_fasta_int_chr_${current_chrom}.bed | bgzip > anc_fasta_intersect/anc_is_ref_chr_${current_chrom}.bed.gz
awk '{if($4!=$5) print $1"\t"$2"\t"$3"\t"$4"\t"$5}' anc_fasta_intersect/anc_fasta_int_chr_${current_chrom}.bed | bgzip > anc_fasta_intersect/anc_not_ref_chr_${current_chrom}.bed.gz

rm anc_fasta_intersect/tmp_anc_fasta_int_chr_${current_chrom}.bed anc_fasta_intersect/anc_fasta_int_chr_${current_chrom}.bed