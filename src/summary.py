import os
import glob
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from sklearn.decomposition import PCA
from sklearn.cross_decomposition import CCA

# Paths
MOD_PATH = "data/embeddings/all_genes/"
INTERSECT_PATH = "data/embeddings/intersect/"
META_FILE = "data/embed_meta.csv"
RESULTS_PATH = "results/plots/"
os.makedirs(RESULTS_PATH, exist_ok=True)


# Collect gene sets from raw embedding data
def collect_gene_sets(mod_path):
    genes_dict = {}
    subdirs = [
        d for d in os.listdir(mod_path) if os.path.isdir(os.path.join(mod_path, d))
    ]

    for subdir in subdirs:
        subdir_path = os.path.join(mod_path, subdir)
        txt_files = [f for f in os.listdir(subdir_path) if f.endswith(".txt")]

        if txt_files:
            txt_path = os.path.join(subdir_path, txt_files[0])
            with open(txt_path, "r") as f:
                genes = set(line.strip() for line in f if line.strip())
                genes_dict[subdir] = genes
    return genes_dict


# Compute Jaccard matrix
def compute_jaccard_matrix(genes_dict):
    subdirs = list(genes_dict.keys())
    jaccard_matrix = np.zeros((len(subdirs), len(subdirs)))

    for i, key1 in enumerate(subdirs):
        for j, key2 in enumerate(subdirs):
            set1 = genes_dict[key1]
            set2 = genes_dict[key2]
            intersection = len(set1 & set2)
            union = len(set1 | set2)
            jaccard_matrix[i, j] = intersection / union if union > 0 else 0

    return pd.DataFrame(jaccard_matrix, index=subdirs, columns=subdirs)


# Plot heatmap
def plot_heatmap(data, title, output_file):
    plt.figure(figsize=(8, 6))
    sns.clustermap(data, annot=False, cmap="viridis", cbar=True, fmt=".2f")
    plt.title(title)
    plt.savefig(output_file, format="pdf", bbox_inches="tight")
    plt.show()


# Process intersect folder for embeddings and gene lists
def process_intersect(intersect_path):
    subfolders = [f.path for f in os.scandir(intersect_path) if f.is_dir()]
    gene_lists = {}
    embeddings = {}

    for subfolder in subfolders:
        gene_txt_files = glob.glob(os.path.join(subfolder, "*.txt"))
        if gene_txt_files:
            with open(gene_txt_files[0], "r") as f:
                genes = [line.strip() for line in f]
            gene_lists[subfolder] = genes

        csv_files = glob.glob(os.path.join(subfolder, "*.csv"))
        if csv_files:
            embeddings[subfolder] = pd.read_csv(csv_files[0], header=None)

    return gene_lists, embeddings


# Align embeddings to the master gene list
def align_embeddings(gene_lists, embeddings):
    master_gene_list = gene_lists[list(gene_lists.keys())[0]]
    embeddings_aligned = {}

    for subfolder, embedding in embeddings.items():
        embedding.index = gene_lists[subfolder]
        embeddings_aligned[subfolder] = embedding.reindex(master_gene_list)

    return embeddings_aligned


# Compute weighted Jaccard similarity
def compute_weighted_jaccard(scaled_correlation_vectors):
    subfolder_names = list(scaled_correlation_vectors.keys())
    n_embeddings = len(subfolder_names)
    weighted_jaccard_matrix = np.zeros((n_embeddings, n_embeddings))

    for i in range(n_embeddings):
        for j in range(i, n_embeddings):
            vec1 = scaled_correlation_vectors[subfolder_names[i]]
            vec2 = scaled_correlation_vectors[subfolder_names[j]]
            min_sum = np.sum(np.minimum(vec1, vec2))
            max_sum = np.sum(np.maximum(vec1, vec2))
            wj_sim = min_sum / max_sum if max_sum != 0 else 0
            weighted_jaccard_matrix[i, j] = wj_sim
            weighted_jaccard_matrix[j, i] = wj_sim

    return pd.DataFrame(
        weighted_jaccard_matrix, index=subfolder_names, columns=subfolder_names
    )


# Clean subfolder names
def clean_subfolder_name(name, base_path="data/embeddings/intersect/", suffix=""):
    if name.startswith(base_path):
        name = name[len(base_path) :]
    if name.endswith(suffix):
        name = name[: -len(suffix)]
    return name

def apply_pca_to_embeddings(embeddings_dict, n_components=10):
    pca_embeddings = {}
    for name, embedding in embeddings_dict.items():
        pca = PCA(n_components=n_components)
        embedding_reduced = pca.fit_transform(embedding)
        pca_embeddings[name] = embedding_reduced
    return pca_embeddings

def compute_similarity(emb1, emb2):
    n_components = min(emb1.shape[1], emb2.shape[1])
    cca = CCA(n_components=n_components)
    cca.fit(emb1, emb2)
    return cca.score(emb1, emb2)

if __name__ == "__main__":
    # Collect gene sets
    genes_dict = collect_gene_sets(MOD_PATH)
    print("Gene sets collected.")

    # Compute Jaccard matrix
    jaccard_df = compute_jaccard_matrix(genes_dict)
    plot_heatmap(
        jaccard_df,
        "Jaccard Index of Genes Heatmap",
        f"{RESULTS_PATH}/jaccard_index_of_genes.pdf",
    )
    print("Jaccard heatmap plotted.")

    # Process intersect folder
    gene_lists, embeddings = process_intersect(INTERSECT_PATH)
    embeddings_aligned = align_embeddings(gene_lists, embeddings)

    # Compute correlation matrices and weighted Jaccard
    gene_correlation_matrices = {
        subfolder: embedding.transpose().corr(method="pearson").fillna(0)
        for subfolder, embedding in embeddings_aligned.items()
    }

    scaled_correlation_vectors = {
        subfolder: (gene_corr.values[np.triu_indices(len(gene_corr), k=1)] + 1) / 2
        for subfolder, gene_corr in gene_correlation_matrices.items()
    }

    wj_sim_df = compute_weighted_jaccard(scaled_correlation_vectors)
    wj_sim_df_clean = wj_sim_df.copy()

    # Clean names and plot
    wj_sim_df_clean.index = wj_sim_df_clean.index.map(clean_subfolder_name)
    wj_sim_df_clean.columns = wj_sim_df_clean.columns.map(clean_subfolder_name)
    plot_heatmap(
        wj_sim_df_clean,
        "Weighted Jaccard Similarity of Gene Correlations",
        f"{RESULTS_PATH}/summary_weighted_heatmap.pdf",
    )
    print("Weighted Jaccard heatmap plotted.")
