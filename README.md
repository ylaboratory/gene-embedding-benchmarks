# Benchmarking gene embeddings from sequence, expression, network, and text models for functional prediction tasks

This repository contains code and datasets for benchmarking gene embedding methods across
individual gene attributes, paired gene interactions, and gene set relationships. This is a
companion resource to our forthcoming study which can be read at the citation below.

>Benchmarking gene embeddings from sequence, expression, network, and text models for functional prediction tasks.
 Zhong J, Li L, Dannenfelser R, and Yao V. bioRxiv (2025)
 [https://doi.org/10.1101/2025.01.29.635607](https://doi.org/10.1101/2025.01.29.635607)

## About

Accurate, data-driven representations of genes are critical for interpreting high-throughput biological data, yet no consensus exists on the most effective embedding strategy for common functional prediction tasks. Here, we present a systematic comparison of 38 gene embedding methods derived from amino acid sequences, gene expression profiles, protein-protein interaction networks, and biomedical literature. We benchmark each approach across three classes of tasks: predicting individual gene attributes, characterizing paired gene interactions, and assessing gene set relationships while trying to control for data leakage. Overall, we find that literature-based embeddings deliver superior performance across prediction tasks, sequence-based models excel in genetic interaction predictions, and expression-derived representations are well-suited for disease-related associations. Interestingly, network embeddings achieve similar performance to literature-based embeddings on most tasks despite using significantly smaller training sets. The type of training data has a greater influence on performance than the specific embedding construction method, with embedding dimensionality having only minimal impact. Our benchmarks clarify the strengths and limitations of current gene embeddings, providing practical guidance for selecting representations for downstream biological applications.

## Organization

This repo is organized into several sections, with the gene embeddings stored on [embedding zenodo](https://zenodo.org/records/16764517).
If rerunning all code from the embedding stage, please download these files and put them under the `data` folder in a 
new folder called `embeddings/intersect` for the common genes and `embeddings/all_genes` for the full embeddings.
- `data`: contains datasets and metadata used for benchmarking
  - `gmt`: gene set files used for benchmarking 
  - `matched_pairs`: files used to map one annotation to another
  - `obo`: ontology files for hierarchical biological relationships
  - `paired_gene_interaction_data`: (downloadable from BioGRID) files used for benchmarking paired genetic interactions
  - `slim_sets`: subsets of annotation terms
  - `embed_meta.csv`: metadata file detailing the embedding methods, their training input type, algorithm, and dimension.
- `results`: houses a folder `tsvs` containing tab delimited output files used to recreate all results in the paper, the `plots` directory is a placeholder directory for figures made after running `paper_figures.R`
- `src`: contains the code used for preprocessing, summarizing, and benchmarking embeddings across our functional prediction tasks
  - `gene_level_benchmark`: code used for benchmarking disease gene prediction (OMIM) and gene function prediction (GO).
  - `gene_pair_benchmark`: code used for benchmarking genetic interaction (e.g., SL/NG) and transcription factor target (TF) prediction.
  - `gene_set_benchmark`: code used for benchmarking matching pathways (GO/KEGG) and disease/tissue (OMIM/Brenda).
  - `paper_figures.R`: code used for reproducing exact figures from the paper (uses the outputs in `results/tsvs`)
  - `other`: code used for preprocessing embeddings and conducting the CCA and ANOVA analysis
    - `preprocess_embedding`: code used for preprocessing embeddings

## Environment and Dependencies

We recommend using conda for installing all necessary packages. Once conda is
installed, get started by creating and activating the virtual environment.

```bash
conda env create -f env.yml
conda activate gene_embed_benchmark 
```

To run the gene set benchmarks, you will also need to install [ANDES](https://pubmed.ncbi.nlm.nih.gov/39231608/). To run ANDES in distinct mode (no overlap), set distinct = True in ANDES/src/set_analysis_func.py in the andes() function.
Please refer to the [ANDES repository](https://github.com/ylaboratory/ANDES) for futher setup and usage instructions. 