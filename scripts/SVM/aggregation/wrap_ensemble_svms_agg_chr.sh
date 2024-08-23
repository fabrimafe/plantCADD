#!/bin/bash
mem=$1
dataset_folder=$2
prediction_folder=$3
name=$(echo "$dataset_folder" | awk -F'_' '{print $1}')

for chr in {1..12}
do
    bsub -n 1 -q new-short -J agg_ensemble_svm_${name}_chr_${chr} -e agg_ensemble_svm_${name}_chr_${chr}_%J.err -o agg_ensemble_svm_${name}_chr_${chr}_%J.out \
    -R rusage[mem=${mem}] "python aggregate_predictions.py df_chr_${chr}_normed.csv df_chr_${chr}_pos.csv ${dataset_folder} ${prediction_folder} ${name}" 
done