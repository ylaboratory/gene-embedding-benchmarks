# Benchmarking gene embeddings from sequence, expression, network, and text models for functional prediction tasks
This repository contains code and datasets for benchmarking gene embedding methods across individual gene attributes, paired gene interactions, and gene set relationships.

## About
Gene embeddings have emerged as transformative tools in computational biology, enabling the efficient translation of complex biological datasets into compact vector representations. This study presents a comprehensive benchmark by evaluating 38 classic and state-of-the-art gene embedding methods across a spectrum of functional prediction tasks. These embeddings, derived from data sources such as amino acid sequences, gene expression profiles, protein-protein interaction networks, and biomedical literature, are assessed for their performance in predicting individual gene attributes, paired gene interactions, and gene set relationships. Our analysis reveals that biomedical literature-based embeddings consistently excel in general predictive tasks, amino acid sequence embeddings outperform in functional and genetic interaction predictions, and gene expression embeddings are particularly well-suited for disease-related tasks. Importantly, we find that the type of training data has a greater influence on performance than the specific embedding construction method, with embedding dimensionality having only minimal impact. By elucidating the strengths and limitations of various gene embeddings, this work provides guidance for selecting and successfully leveraging gene embeddings for downstream biological prediction tasks.

## Organization
This repo is organized into several sections, part of which is stored on [zenodo](https://zenodo.org/records/14769058).
- `bin`: contains binaries and intermediate files from the benchmarking experiments, which includes the fold and holdout splits that we used in our tests saved as pkl files
- `data`: contains datasets and metadata used for benchmarking
  - `embeddings`: preprocessed embeddings for genes from various methods
    - `intersect`: preprocessed embeddings for genes that are common across all methods in entrez gene format (on [zenodo](https://zenodo.org/records/14769058))
    - `all_genes`: preprocessed embeddings that contain all genes in entrez gene format
  - `gmt`: gene set files used for benchmarking 
  - `matched_pairs`: files used to map one annotation to another
  - `obo`: ontology files for hierarchical biological relationships
  - `paired_gene_interaction_data`: files used for benchmarking paired genetic interactions (downloadable from BioGRID)
  - `slim_sets`: subsets of annotation terms
  - `embed_meta.csv`: metadata file detailing the embedding methods, their training input type, algorithm, and dimension.
- `results`: contains the results of the gene level and gene pair benchmarking experiments
  - `andes_results`: contains the scores from the gene set benchmarks (on [zenodo](https://zenodo.org/records/14769058))
- `src`: contains the code used for preprocessing, summarizing, and benchmarking embeddings across our functional prediction tasks
  - `gene_level_benchmark`: code used for benchmarking disease gene prediction (OMIM) and gene function prediction (GO).
  - `gene_pair_benchmark`: code used for benchmarking genetic interaction (e.g., SL/NG) and transcription factor target (TF) prediction.
  - `gene_set_benchmark`: code used for benchmarking matching pathways (GO/KEGG) and disease/tissue (OMIM/Brenda).
  - `preprocess_embedding`: code used for preprocessing embeddings
  - `summary.py`: code used for summarizing the tested embeddings


## Running the scripts
We recommend using conda for installing all necessary packages. Once conda is installed, get started by creating and activating the virtual environment.

 ```bash
 conda env create -f env.yml
 conda activate gene_embed_benchmark 
 ```
 To run the gene set benchmarks, first download the ANDES tool:
```bash 
git clone https://github.com/ylaboratory/ANDES.git
```
Make sure to navigate to the appropriate directory and follow any additional instructions provided in the [ANDES repository](https://github.com/ylaboratory/ANDES) for setting up and running the tool.

Each folder in the `src` directory contains either a Python (.py) or Jupyter Notebook (.ipynb) file that corresponds to a specific benchmark or analysis. These scripts not only run the benchmarks but also generate plots, which are saved in the `results` folder. To execute the scripts, use the following commands:

For Python files:
```bash
python path_to_script/script_name.py
```
For Jupyter Notebooks, open them in a Jupyter environment:
```bash
jupyter notebook path_to_script/script_name.ipynb
```
Most of the scripts rely on their corresponding helper.py file, which contains helper functions used throughout the analyses. 

## Citation
>Benchmarking gene embeddings from sequence, expression, network, and text models for functional prediction tasks.
 Zhong J, Li L, Dannenfelser R, and Yao V. bioRxiv (2025)
 [https://doi.org/10.1101/2025.01.29.635607](https://doi.org/10.1101/2025.01.29.635607)
