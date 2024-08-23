#!/bin/bash
datasets_folder=$1
prediction_folder=$2
mem=$3
for i in {0..49}
do
    bsub -n 1 -q new-short -J predict_ensemble_svm_${i}_test -e predict_ensemble_svm_${i}_test_%J.err -o predict_ensemble_svm_${i}_test_%J.out \
    -R rusage[mem=$mem] "python predict_mini_SVM.py ${i} all_4.28M_test_df_normed.csv ${datasets_folder} ${prediction_folder}" 
done