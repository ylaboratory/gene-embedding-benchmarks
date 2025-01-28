# Benchmarking gene embeddings from sequence, expression, network, and text models for functional prediction tasks
This repository contains code and datasets for benchmarking 38 gene embedding methods across functional prediction tasks, including gene attribute, interaction, and set relationship predictions, using data from amino acid sequences, gene expression, PPI networks, and biomedical literature.

## Organization
This repo is organized into several sections, part of which is stored on zenodo.
- `data`: contains datasets and metadata used for benchmarking
- `src`: contains the code used for related tasks to preprocess the embeddings and benchmark them across our three main functional prediction tasks
The data section is further divided into sections:
- `embeddings`: preprocessed embeddings for genes from various methods (on zenodo)
- `gmt`: gene set files used for benchmarking 
- `matched_pairs`: files used to map one annotation to another
- `obo`: ontology files for hierarchical biological relationships
- `slim_sets`: subsets of annotation terms
- `embed_meta.csv`: metadata file detailing the embedding methods, their training input type, algorithm, and dimension.

## Running the scripts
We recommend using conda for installing all necessary packages. Once conda is installed 
get started by creating and activating the virtual environment.

 ```bash
 conda env create -f env.yml
 conda activate gene_embed_benchmark 
 ```
 To run the gene set benchmarks, first clone the necessary repository:
```bash 
git clone https://github.com/ylaboratory/ANDES.git
```
Make sure to navigate to the appropriate directory and follow any additional instructions provided in the [ANDES repository](https://github.com/ylaboratory/ANDES) for setting up and running the tool.

## Citation
