#!/bin/bash
# Modules
module load tabix
module load miniconda/22.11.1_environmentally
conda activate NN_new

# Args
gff_file=$1
bed_file=$2
chromosome=$3

filename=$(basename "$gff_file")
org_name="${filename%%.*}"

# Get distance
python calc_gene_dist.py ${gff_file} ${bed_file} ${chromosome} | awk '{ gsub("inf", "-1"); print }' | bgzip > ${org_name}_both_gene_dist_chr_${chromosome}.bed.gz