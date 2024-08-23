#!/bin/bash

nucleotide_file_name=$1
chr=$2
cat fasta_nucleotides/${nucleotide_file_name}_chr${chr}.bed |awk '{ $4 = toupper($4); print $1"\t"$2"\t"$3"\t"$4}' | sort -k1,1n -k2,2n > fasta_nucleotides/${nucleotide_file_name}_parsed_chr${chr}.bed
rm fasta_nucleotides/${nucleotide_file_name}_chr${chr}.bed