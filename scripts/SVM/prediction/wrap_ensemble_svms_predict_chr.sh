#!/bin/bash

chr=$1
mem=$2
datasets_folder=$3
prediction_folder=$4
name=$(echo "$datasets_folder" | awk -F'_' '{print $1}')

for i in {0..49}
do
    bsub -n 1 -q new-short -J predict_ensemble_svm_${i}_${name}_chr_${chr} -e predict_ensemble_svm_${i}_${name}_chr_${chr}_%J.err -o predict_ensemble_svm_${i}_${name}_chr_${chr}_%J.out \
    -R rusage[mem=$mem] "python predict_mini_SVM.py ${i} df_chr_${chr}_normed.csv ${datasets_folder} ${prediction_folder}" 
done