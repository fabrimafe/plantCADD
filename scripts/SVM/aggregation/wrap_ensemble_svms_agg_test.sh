#!/bin/bash

bsub -n 1 -q new-short -J agg_ensemble_svm_test -e agg_ensemble_svm_test_%J.err -o agg_ensemble_svm_test_%J.out \
-R rusage[mem=3000] "python aggregate_predictions.py all_4.28M_test_df_normed.csv 4.28M_test_df_positions.csv test_datasets test_predictions"