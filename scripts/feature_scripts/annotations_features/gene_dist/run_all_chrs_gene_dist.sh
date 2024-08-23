#!/bin/bash
# A script to run all gene dist chromosomes at once

#chromosomes=$1 # The main script will se#nd this
gff_file=$1
bed_file_prefix=$2
memory=$3

for i in {0..12} # The main script will send this, for now, only to make sure it works, manually.
do
    bsub -n 1 -m public_hosts -q new-long -J generate_gene_dist_chr_${i} -e generate_gene_dist_chr_${i}.error -o generate_gene_dist_chr_${i}.output -R rusage[mem=${memory}] "bash generate_gene_dist.sh ${gff_file} ${bed_file_prefix}_chr_${i}.bed.gz ${i}"
done