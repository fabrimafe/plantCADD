#!/bin/bash
# This script is a wrapper to check the enrichment for the neutral mutation dataset using the ancestral state files fabrizio creates
# The reference genome should be in the general folder
# The ancestral file should also be in the general folder, gzipped
# This script should be in a scripts folder inside the general folder

# There are three groups to check
# 1. Completely fixed, where ancestral != reference
# For this one we need a file where we have ancestral and reference states
# For that we need reference states
ref_genome=$1
anc_file=$2
chrom_list=$3


# Get the fasta in bed format
echo "Parsing fasta to bed..."
# parameters are the fasta file and output directory
python scripts/parse_fasta_to_bed_chroms.py ${ref_genome} fasta_nucleotides

# Also split the anc file into chromosomes for later
cd anc_states
echo "Splitting ancestral state file to chromomsomes..."
zcat "../${anc_file}" | awk '{print > $1}'

# Change the chromosome names into bed format
for i in ${chrom_list}
do
  mv ${i} anc_state_chr_${i}.bed
done

# Back to main folder
cd ..