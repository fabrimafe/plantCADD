#!/bin/bash
# The main bash script to run all the other subscripts
# THE FIRST SCRIPT ASSUMES THE FASTA AND ANC (BGZIPPED) FILES ARE IN THE MAIN FOLDER
# THE SECOND SCRIPT ASSUMES THE ANC STATE HAS 4 COLS AND SO DOES THE FASTA NUC
# THE THIRD SCRIPT ASSUMES THE VCF FILES ARE IN THE ./vcfs FOLDER. THIER NAMES SHOULD BE 'CHR_NAME.vcf.gz'
# THE FOURTH SCRIPT ASSUMES A GFF FILE WITH CDS POSITIONS IN THE MAIN FOLDER

# For quiet bsub commands for this script

# Modules
module load bedtools
conda activate NN_new

# Parameters
ref_genome=$1
anc_state=$2
gff_file=$3
home_dir=$(pwd)

# Making sure subdirectories exist
mkdir -p fasta_nucleotides
mkdir -p anc_states
mkdir -p anc_fasta_intersect
mkdir -p completely_fixed
mkdir -p nearly_fixed
mkdir -p rare_mutations
mkdir -p codon_pos
mkdir -p enrichment

# ----------------------------------PARSING OF FILES------------------------------

# VCF parsing
echo "Parsing VCF files..."
for i in ${chromosomes}
do
    bsub -n 1 -m public_hosts -q new-short -J ${i}_parse_vcf -R rusage[mem=500] "bash scripts/3_parse_vcf.sh ${i}"
done

# Creating codon position files
# Generate the specific codon position
echo "Parsing codon position..."
bash scripts/4_codon_pos.sh ${gff_file} ${chromosomes}