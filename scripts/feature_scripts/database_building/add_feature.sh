#!/bin/bash
# A script to add a feature to the feature file
# To add a few columns, change feature_num
# When appending features to the feature file, from the main script, I need to make sure there are no headers in the files
# This is for GERP and PHYLOP specifically, my own features are header-less.
export TMPDIR=/scratch/

module load bedtools #; module load BEDTools
module load tabix #;module load tabixpp

chr=$1
feature_name=$2
org_name=$3
map_or_int=$4
column_to_add=$5
missing_value=$6

feature_file=/home/labs/alevy/omerbar/features/arabidopsis/main_table/${org_name}_features_SIFT_chr_Chr${chr}.bed.gz
new_feature=/home/labs/alevy/omerbar/expanded_alignment_sl5/PRIVATE_WEXAC_QUEUES_RAXML/ARAB_ALL/results/gerp/Chr${chr}.gerp.final.bed.gz
feature_file_lines=$(wc -l < <(zcat $feature_file))

echo $map_or_int
if [ "${map_or_int}" = "map" ]; then
  echo "Mapping..."
  bedtools map -a <(zcat $feature_file) -b <(zcat $new_feature) -c ${column_to_add} -null ${missing_value}  > tmp_intersected_result_chr_${chr}.bed # | sed 's/\r//g'
elif [ "${map_or_int}" = "intersect" ]; then
  echo "Intersecting..."
  bedtools intersect -a <(zcat $feature_file) -b <(zcat $new_feature) -sorted -wa -wb > tmp_intersected_result_chr_${chr}.bed # | sed 's/\r//g'
  feature_file_cols=$(zcat $feature_file | awk 'NR==1{print NF; exit}')
  awk -v num_cols="$feature_file_cols" \
  -v OFS="\t" '{
    # Print all columns from intersected_result.bed 
    for (i = 1; i <= num_cols; i++) printf "%s\t", $i

    # Print the last X columns from the second BED file
    for (i = num_cols + 4; i <= NF; i++) {
      printf "%s", $i
      if (i < NF) printf "\t"
    }

    printf "\n"
  }' tmp_intersected_result_chr_${chr}.bed > tmp_clean_result_chr_${chr}.bed
  mv tmp_clean_result_chr_${chr}.bed tmp_intersected_result_chr_${chr}.bed
else
  echo "Choose map or intersect."
  exit
fi

echo "Validating..."
# Check the lines of this file
new_file_lines=$(wc -l < "tmp_intersected_result_chr_${chr}.bed")
# If we lost lines, i.e. files aren't the same
if [ "$new_file_lines" -ne "$feature_file_lines" ]; then
  echo "Something went wrong, number of lines not the same. Exiting..."
  echo "Feature file: ${feature_file_lines}"
  echo "New file: ${new_file_lines}"
  exit
else
  echo "Wipply doo!"
  bgzip "tmp_intersected_result_chr_${chr}.bed" 
  mv "tmp_intersected_result_chr_${chr}.bed.gz" $feature_file
fi
