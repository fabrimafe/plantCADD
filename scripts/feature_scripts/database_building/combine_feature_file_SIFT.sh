#!/bin/bash
# A script to merge the new SIFT scores with the feature file.
# The structure of the feature file is changing, allowing for multiple rows for the same nucleotide for CDS regions.
# This will be acheived in two parts, one of which will be to get all the positions for which we have mutations,
# i.e. CDS regions. This will be acheieved with bedtools intersect.
# The second part would be to take all the other nucleotides and append them, since those do not have any sort of information over them (SIFT info), we will simply add 
# them and append the information. 
# For wexac jobs
export TMPDIR=/scratch/
# Modules
module load bedtools
module load tabix

# Args
chr=$1
org_name=$2

feature_file=/home/labs/alevy/omerbar/features/arabidopsis/main_table/${org_name}_features_chr_${chr}.bed.gz
SIFT_file=/home/labs/alevy/omerbar/features/arabidopsis/SIFT/annotations/FINAL_SIFT_predictions_${chr}.bed.gz

# Get number of cols of feature file
num_cols_feat=$(awk '{print NF; exit}' <(zcat ${feature_file}))

# Step 1
bedtools intersect -a <(zcat ${feature_file}) -b <(zcat ${SIFT_file}) -wa -wb | \
awk -v num_cols="$num_cols_feat" '{ for (i = 1; i <= num_cols; i++) printf "%s\t", $i; for (i = num_cols + 5; i < NF; i++) printf "%s\t", $i; printf "%s\n", $NF }'  >  tmp_${org_name}_chr_${chr}_CDS.bed
 
# Step 2
bedtools subtract -a <(zcat ${feature_file}) -b  tmp_${org_name}_chr_${chr}_CDS.bed | \
awk -v num_cols="$num_cols_feat" '{ for (i = 1; i <= num_cols; i++) printf "%s\t", $i; for (i = num_cols + 1; i <= num_cols + 4; i++) { printf ".\t"}; for (i = 0; i <= 3; i++) printf "-1\t"; printf "%s\n", "TOLERATED" }' > tmp_${org_name}_chr_${chr}_NON_CDS.bed

# Create a header file
#echo '#CHR\tSTART\tEND\tREF_ALLELE\tKMERS_20\tGC_35\tCODON_POS\tGEN_REGION\tDIST_UP\tDIST_DOWN\tGERP_ExpSubst\tGERP_RejSubstScore\tGERP_TaxaAligned\tPHYLOP\tALT_ALLELE\tREF_AMINO\tALT_AMINO\tVARIANT_TYPE\tAMINO_POS\tSIFT_SCORE\tSIFT_MEDIAN\tNUM_SEQS\tSIFT_PREDICTION' > /home/labs/alevy/omerbar/features/expanded_SL5/main_table/${org_name}_features_SIFT_chr_${chr}.bed
#
# Concatenate them together and sort
cat tmp_${org_name}_chr_${chr}_CDS.bed tmp_${org_name}_chr_${chr}_NON_CDS.bed | \
sort -k1,1n -k2,2n >> /home/labs/alevy/omerbar/features/arabidopsis/main_table/${org_name}_features_SIFT_chr_${chr}.bed

# Bgzip
bgzip /home/labs/alevy/omerbar/features/arabidopsis/main_table/${org_name}_features_SIFT_chr_${chr}.bed

# Remove tmps
#rm tmp_${org_name}_chr_${chr}_NON_CDS.bed tmp_${org_name}_chr_${chr}_CDS.bed