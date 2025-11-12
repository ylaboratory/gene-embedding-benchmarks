#!/bin/bash

# Input arguments
SUBFOLDER="path/to/embedding_subfolder"  
OUT_ROOT="path/to/output_directory"

CV_FOLD1_PKL="data/data_splits/gene_level_benchmark/go_folds_splits/go_cv_fold1_dict_all.pkl"
CV_FOLD2_PKL="data/data_splits/gene_level_benchmark/go_folds_splits/go_cv_fold2_dict_all.pkl"
CV_FOLD3_PKL="data/data_splits/gene_level_benchmark/go_folds_splits/go_cv_fold3_dict_all.pkl"
HOLDOUT_PKL="data/data_splits/gene_level_benchmark/go_folds_splits/go_holdout_dict_all.pkl"

python gene_level_benchmarks.py \
  --subfolder "$SUBFOLDER" \
  --cv-fold1-pkl "$CV_FOLD1_PKL" \
  --cv-fold2-pkl "$CV_FOLD2_PKL" \
  --cv-fold3-pkl "$CV_FOLD3_PKL" \
  --holdout-pkl "$HOLDOUT_PKL" \
  -d "$OUT_ROOT"
