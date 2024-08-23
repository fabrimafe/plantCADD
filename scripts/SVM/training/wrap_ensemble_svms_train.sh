#!/bin/bash

mem=$1
batch_size=$2
desc=$3

for classifier_id in {0..49}
do
    bsub -n 1 -q new-short -J train_ensemble_svm_${desc}_${classifier_id} -e train_ensemble_svm_${desc}_${classifier_id}_%J.err -o train_ensemble_svm_${desc}_${classifier_id}_%J.out\
    -R rusage[mem=${mem}] "python train_mini_SVM.py bootstrap_sample_${classifier_id}.csv.gz ${batch_size} ${classifier_id}" 
done
