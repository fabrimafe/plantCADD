

#!/bin/bash
name_file=$1
mem=$2
target_folder=$3
mkdir -p ${target_folder}/scripts
mkdir -p ${target_folder}/logs

# Go over genome names from genome file, given as a parameter to this script
while IFS= read -r line; do
  # Ceeate a script for the bjob
  cd ${target_folder}/scripts
  echo "##!/bin/bash
# This file is used to send the msa_pipeline work into the cluster.
# To do this I need to load some modules:
module load miniconda/22.11.1_environmentally
module load snakemake

# I initialize the conda workspace
#conda init bash
#conda init bash
# And activate the anaconda env 'msa':
conda activate msa

# Next I navigate into the pipeline's folder
cd ${target_folder}/${line}_pipeline
# And run the first pipeline command:
snakemake --use-conda --rerun-incomplete -j 99 --configfile config/config.yaml -R align" > "${line}_align.sh"

  # Run the bjob from the logs folder, to have logs there  
  cd ${target_folder}/logs
  bsub -n 1 -q new-long -m public_hosts -J  pairwise_align_${line} -e  pairwise_align_${line}_%J.error -o pairwise_align_${line}_%J.output -R rusage[mem=${mem}] ${target_folder}/scripts/${line}_align.sh 
  #sleep 3
done < ${name_file}