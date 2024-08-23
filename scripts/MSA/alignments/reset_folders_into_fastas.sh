

#!/bin/bash
name_file=$1
target_folder=$2
mkdir -p ${target_folder}/genome_fastas
# Go over genome names from genome file, given as a parameter to this script
while IFS= read -r line; do
  echo "${line}"
  # Ceeate a script for the bjob
  mv ${target_folder}/${line}_pipeline/data/${line}.fa ${target_folder}/genome_fastas/${line}.fa
  rm -rf ${target_folder}/${line}_pipeline &
done < ${name_file}