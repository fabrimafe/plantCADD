#!/bin/bash
name_file=$1
target_folder=$2
# Go over genome names from genome file, given as a parameter to this script
while IFS= read -r line; do
  # Just debug
  echo "${line}"
  cd ${target_folder}
  # Create the pipeline directory for this genome
  mkdir "./${line}_pipeline"
  # Go into the directory to clone the git
  cd "${line}_pipeline"
  # Clone the git
  git clone https://bitbucket.org/bucklerlab/msa_pipeline .
  # Add config file
  echo "##############################################
# Example config for msa_pipeline Snakefile  #
##############################################

########################
# REQUIRED PARAMETERS  #
########################

# Reference species name
# The name should match the name of a fasta file in "msa_pipeline/data" but NOT contain the .fa suffix
refName: A_thaliana

# Query species to be aligned to reference
# The species names should each match a fasta file in "msa_pipeline/data" but NOT contain the .fa suffix
species:
  - ${line}

######################################
# OPTIONAL PARAMETERS  FOR ALIGNMENT #
######################################

# Split input fasta files into N chunks for improved parallel processing
# Chunks are aligned in parallel if sufficient threads are available
# Setting N to be greater than the minimum chromosome/scaffold number may lead to errors
splitFastaN: 1

# Set alignment tool: last|minimap2|gsalign
aligner: last

# Change default alignment parameters
lastParams: \"-m 10 -j 3 -u 1 -p HOXD70\"
minimap2Params: \"-a -cx asm20\"
gsalignParams: \"-sen -no_vcf\"

# Without last-split, the alignments are many-many after lastal, and many-one afer
# chaining and netting.
lastSplit:
# For one-to-one last alignments, comment the above and uncomment the lastSplit line below
#lastSplit: \" | last-split | maf-swap | last-split | maf-swap \"

# Roast multiple alignment parameters
roastParams: \"+ X=2 E=\"

# Newick format species tree for roast
# It must contain the reference and all other species in the alignment
# If left empty, tree will be calculated from genome sequences with mashtree
speciesTree: 

#################################################
# OPTIONAL PARAMETERS FOR CONSERVATION ANALYSIS #
#################################################

# When callling conservation set the max number of sites to use to calculate the neutral model
# Higher numbers may slow down the analysis without providing a better model
maxNeutralSites: 10000000

# A GFF file with CDS features must be provided in the /data dir to use a neutral model based on 4-fold degenerate sites
refGFF: S_lycopersicum_heinz.SL4" > ./config/config.yaml
  # Go back to main folder(~/all_genomes/genomes)
#  cd ../..
  # Move the corresponding fasta file for this genome into the pipeline folder
  cp "${target_folder}/genome_fastas/${line}.fa" "./data/${line}.fa" &
  cp /home/labs/alevy/omerbar/backups/A_thaliana.fa "./data/A_thaliana.fa" &
  # Copy the referfence genome (Lycopersicum) into the pipeline data folder
#  cp S_lycopersicum_heinz.fa "${line}_pipeline/data"
  # Remove executable just in case
 # chmod -R a-x "${line}_pipeline/data/S_lycopersicum_heinz.fa"
  # Move the folder to the pipeline folders
  #mv "${line}_pipeline" pipeline_folders
done < ${name_file}