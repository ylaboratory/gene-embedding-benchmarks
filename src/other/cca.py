from sklearn.decomposition import PCA
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cross_decomposition import CCA
import matplotlib.patches as mpatches
import pandas as pd
import glob
import os


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
    # X_c, Y_c = cca.transform(emb1, emb2)
    return cca.score(emb1, emb2)


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


MOD_PATH = "data/embeddings/all_genes/"
INTERSECT_PATH = "data/embeddings/intersect/"

gene_lists, embeddings = process_intersect(INTERSECT_PATH)
embeddings_aligned = align_embeddings(gene_lists, embeddings)

embeddings_aligned_np = {
    name: np.array(embedding) for name, embedding in embeddings_aligned.items()
}
embeddings_aligned_np_pca = apply_pca_to_embeddings(
    embeddings_aligned_np, n_components=10
)

embedding_names = sorted(embeddings_aligned_np_pca.keys())
n_embeddings = len(embedding_names)
similarity_matrix = np.zeros((n_embeddings, n_embeddings))

for i, name_i in enumerate(embedding_names):
    print(name_i)
    for j, name_j in enumerate(embedding_names[i:], start=i):
        sim = compute_similarity(
            embeddings_aligned_np_pca[name_i], embeddings_aligned_np_pca[name_j]
        )
        similarity_matrix[i, j] = sim
        similarity_matrix[j, i] = sim

sim_df = pd.DataFrame(similarity_matrix, index=embedding_names, columns=embedding_names)


g = sns.clustermap(
    sim_df,
    annot=False,
    cmap="vlag",
    fmt=".4f",
    square=True,
    linewidths=0,
    cbar_kws={"shrink": 0.5},
    figsize=(20, 20),
)

g.fig.suptitle("Canonical Correlation Analysis Scores", fontsize=16)

plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90, ha="right", fontsize=20)
plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize=20)


plt.savefig("y_summary_cc_weighted_heatmap.pdf", format="pdf", bbox_inches="tight")

plt.show()
