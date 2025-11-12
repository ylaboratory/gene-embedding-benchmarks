#!/bin/bash

# Input arguments
SUBFOLDER="path/to/embedding_subfolder"       
OUT_ROOT="path/to/output_directory"
OPERATION="sum" # options: sum, product, concat
SUFFIX="sum_intersected"   

CV_PKL="data/data_splits/gene_pair_benchmark/ng_nested_cv_splits.pkl" 

python gene_pair_benchmarks.py \
  --subfolder "$SUBFOLDER" \
  --cv-pkl "$CV_PKL" \
  --operation "$OPERATION" \
  --out-root "$OUT_ROOT" \
  --suffix "$SUFFIX"
