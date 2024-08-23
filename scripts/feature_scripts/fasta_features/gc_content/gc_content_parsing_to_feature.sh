#!/bin/bash
# A script to get the gc content per position in the genome
# This script requires the fasta divided into chromosomal files
# Haven't used this yet, maybe needs to be tweaked.
fasta=$1
window_size=$2 # Must be an odd number
nucleotide_file_location=$3 # The folder in which the bed files containing the nucleotide identity per position.
nuclotide_file_prefix=$4 # The perfix for these files. Say SL5_chr_1.bed, write SL5

# If it is not odd
if (($window_size % 2 == 0)); then
echo "Window size must be odd. Terminating."
exit 1
else # If it is odd
# Modules
module load bedtools
module load tabix
module load miniconda/22.11.1_environmentally
conda activate NN_new

# Getting chromosome sizes
chromosomes=$( cat ${fasta} | grep -oE '^>[^[:space:]]+' | sed 's/^>//' )

# Getting chromosome sizes
echo "Getting chromosome sizes..."
python get_chromosome_sizes.py ${fasta}

# Splitting fasta into chromosomal fasta files
mkdir -p fasta_chromosomal_files
cd fasta_chromosomal_files
awk '/^>/{filename=sprintf("%s.fa", substr($0,2));}{print > filename;}' ${fasta}
cd ..

half_window_size=$(($window_size / 2))
for i in ${chromosomes} # per chromosome
do
	echo "Chromosome ${i}..."
	echo "Making the windows bed file..."
	# Create a window file. This file creates a bed file with ranges per row, describing a 35 window ( when possible ) around any position in the genome. 
	bedtools slop -l ${half_window_size} -r ${half_window_size} -s -i <(zcat ${nucleotide_file_location}/${nuclotide_file_prefix}_chr_${i}_CAP.bed.gz) -g ./chromosomes.sizes > ./${nuclotide_file_prefix}_windows_chr_${i}.bed

	echo "Making the nucleotide bed file..."
	# Make the nucleotide count per window from before, and add it to a file.
    bedtools nuc -fi fasta_chromosomal_files/${i}.fa -bed ./${nuclotide_file_prefix}_windows_chr_${i}.bed > ./${nuclotide_file_prefix}_nuc_content_chr_${i}.bed	
	
	echo "Expanding the bed file..."
	# Make bed format, remove all spaces, zip
	awk 'BEGIN{n=0}NR>1{print $1,"\t",n,"\t",n+1,"\t",$5; n++}' ./${nuclotide_file_prefix}_nuc_content_chr_${i}.bed | sed 's/ //g' | bgzip > ./${nuclotide_file_prefix}_gc_content_window_${window_size}_chr_${i}.bed	 # skip header line with NR>1

	# Remove tmps
	rm ./${nuclotide_file_prefix}_windows_chr_${i}.bed ./${nuclotide_file_prefix}_nuc_content_chr_${i}.bed
done
fi
