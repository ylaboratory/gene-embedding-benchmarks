# Benchmarking gene embeddings from sequence, expression, network, and text models for functional prediction tasks

This repository contains code and datasets for benchmarking gene embedding methods across
individual gene attributes, paired gene interactions, and gene set relationships. This is a
companion resource to our forthcoming study which can be read at the citation below.

>Benchmarking gene embeddings from sequence, expression, network, and text models for functional prediction tasks.
 Zhong J, Li L, Dannenfelser R, and Yao V. bioRxiv (2025)
 [https://doi.org/10.1101/2025.01.29.635607](https://doi.org/10.1101/2025.01.29.635607)

## About

Gene embeddings have emerged as transformative tools in computational biology, enabling the efficient
translation of complex biological datasets into compact vector representations. This study presents a
comprehensive benchmark by evaluating 38 classic and state-of-the-art gene embedding methods across a
spectrum of functional prediction tasks. These embeddings, derived from data sources such as amino acid
sequences, gene expression profiles, protein-protein interaction networks, and biomedical literature,
are assessed for their performance in predicting individual gene attributes, paired gene interactions,
and gene set relationships. Our analysis reveals that biomedical literature-based embeddings consistently
excel in general predictive tasks, amino acid sequence embeddings outperform in functional and genetic
interaction predictions, and gene expression embeddings are particularly well-suited for disease-related
tasks. Importantly, we find that the type of training data has a greater influence on performance than
the specific embedding construction method, with embedding dimensionality having only minimal impact.
By elucidating the strengths and limitations of various gene embeddings, this work provides guidance
for selecting and successfully leveraging gene embeddings for downstream biological prediction tasks.

## Organization

This repo is organized into several sections, with the gene embeddings stored on [embedding zenodo](https://zenodo.org/records/16764517).
If rerunning all code from the embedding stage, please download these files and put them under the `data` folder in a 
new folder called `embeddings/intersect` for the common genes and `embeddings/all_genes` for the full embeddings.
- `bin`: pkl files, which includes the fold and holdout splits that we used in our tests
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
  - `preprocess_embedding`: code used for preprocessing embeddings
  - `summary.py`: code used for summarizing the tested embeddings
  - `paper_figures.R`: code used for reproducing exact figures from the paper (uses the outputs in `results/tsvs`)


## Environment and Dependencies

We recommend using conda for installing all necessary packages. Once conda is
installed, get started by creating and activating the virtual environment.

```bash
conda env create -f env.yml
conda activate gene_embed_benchmark 
```

To run the gene set benchmarks, you will also need to install [ANDES](https://pubmed.ncbi.nlm.nih.gov/39231608/).
Please refer to the [ANDES repository](https://github.com/ylaboratory/ANDES) for futher setup and usage instructions.