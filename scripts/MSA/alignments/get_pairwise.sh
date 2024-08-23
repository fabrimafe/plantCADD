

#!/bin/bash
name_file=$1
target_folder=$2
mkdir -p ${target_folder}/pairwise_alignments


# Go over genome names from genome file, given as a parameter to this script
while IFS= read -r line; do
  # Ceeate a script for the bjob
  mv ${target_folder}/${line}_pipeline/results/toast/A_thaliana.${line}.toast2.maf ${target_folder}/pairwise_alignments/
done < ${name_file}