#!/bin/bash
# Get all fastas from selected pipeline folders

name_file=$1
target_folder=$2
mkdir -p ${target_folder}/ALL/results
mkdir -p ${target_folder}/ALL/results/roast
cp /home/labs/alevy/omerbar/backups/A_thaliana.fa ${target_folder}/ALL/data &

while IFS= read -r line; do
    touch ${target_folder}/ALL/data/${line}.fa
    cp pairwise_alignments/A_thaliana.${line}.toast2.maf ${target_folder}/ALL/results/roast/ &
done < ${name_file}
