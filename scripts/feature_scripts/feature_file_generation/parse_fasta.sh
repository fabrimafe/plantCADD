#!/bin/bash
# A script to extract sequences from the fasta file and put in bed format, divided into chromosomes.
module load tabix
module load miniconda/22.11.1_environmentally
conda activate NN_new

fasta=$1
org_name=${fasta%.*} # Get organism name without .fa prefix
chromosomes=$( cat ${fasta} | grep -oE '^>[^[:space:]]+' | sed 's/^>//' ) # Get chromosome names

# Get the bed file from the fasta file
python parse_fasta_to_bed.py ${fasta} "tmp_${org_name}.bed"

# Sort the bed file
cat "tmp_${org_name}.bed"| \
awk 'BEGIN{OFS="\t"}{$4 = toupper($4)} print' | \ 
sort -k1,1n -k2,2n  > "${org_name}.bed"

# Split it by chromosome
awk '{print > $1}' "${org_name}.bed"

# Change the names of the chromosomes
for i in ${chromosomes}
do
    mv ${i} ${org_name}_chr_${i}.bed
    bgzip ${org_name}_chr_${i}.bed
done