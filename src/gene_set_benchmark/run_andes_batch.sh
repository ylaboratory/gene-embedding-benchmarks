#!/bin/bash

# Input arguments
EMBED_PATH="data/embeddings/intersect" # folder of folders with all embeddings to be tested
ANDES_SCRIPT="ANDES/src/andes.py"
# please install ANDES via the directions in the README
# to run ANDES in distinct mode (no overlap), set distinct = True in ANDES/src/set_analysis_func.py in the andes() function.
GENESET1="data/gmt/KEGG_CPDB.gmt"
GENESET2="data/gmt/hsa_low_eval_BP_propagated.gmt"
OUT_DIR="results/andes_out"
STATUS_LOG="andes_status.txt"
MATCHED_PAIRS="data/matched_pairs/matched_pair_kegg_go.txt"

python batch_andes.py \
  --embed-path "$EMBED_PATH" \
  --andes-script "$ANDES_SCRIPT" \
  --geneset1 "$GENESET1" \
  --geneset2 "$GENESET2" \
  --out-dir "$OUT_DIR" \
  --status-log "$STATUS_LOG" \
  --matched-pairs "$MATCHED_PAIRS"
