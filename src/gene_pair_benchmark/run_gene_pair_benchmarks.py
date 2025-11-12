import pandas as pd
import helper
import pickle
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.api as sm
import statsmodels.formula.api as smf


if __name__ == "__main__":
    embeddings, reference_node2index = helper.load_embeddings()

    # SL
    data = pd.read_csv(
        "/data/biogrid/BIOGRID-ORGANISM-Homo_sapiens-4.4.240.tab3.txt", sep="\t"
    )
    data = data[
        (data["Organism ID Interactor A"] == 9606)
        & (data["Organism ID Interactor B"] == 9606)
    ]
    data = data[data["Experimental System Type"] == "genetic"]

    selected_data = data[data["Experimental System"] == "Synthetic Lethality"]
    positive_pairs = []
    for i in range(selected_data.shape[0]):
        a = str(selected_data.iloc[i]["Entrez Gene Interactor A"])
        b = str(selected_data.iloc[i]["Entrez Gene Interactor B"])
        positive_pairs.append((min(a, b), max(a, b)))
    positive_pairs = set(positive_pairs)

    positive_pairs = {
        pair
        for pair in positive_pairs
        if pair[0] in reference_node2index and pair[1] in reference_node2index
    }
    print(len(positive_pairs))
    used_nodes = list(
        set([x for x, y in positive_pairs]).union(set([y for x, y in positive_pairs]))
    )

    holdout_nodes, splits, cv_nodes = helper.fold_split(used_nodes)
    (
        fold_splits,
        holdout_pairs,
        holdout_labels,
        final_train_pairs,
        final_train_labels,
    ) = helper.setup_data(positive_pairs, splits, holdout_nodes, cv_nodes)

    # save holdout and fold splits
    file_names = [
        "bin/sl_fold_splits_pairs.pkl",
        "bin/sl_fold_nodes.pkl",
        "bin/sl_holdout_pairs.pkl",
        "bin/sl_holdout_labels.pkl",
        "bin/sl_holdout_nodes.pkl",
    ]
    data_dicts = [fold_splits, splits, holdout_pairs, holdout_labels, holdout_nodes]

    for file_name, data_dict in zip(file_names, data_dicts):
        with open(file_name, "wb") as f:
            pickle.dump(data_dict, f)

    sl_fold_results_df, sl_holdout_results_df = helper.run_SVM(
        embeddings,
        reference_node2index,
        splits,
        fold_splits,
        holdout_pairs,
        holdout_labels,
        final_train_pairs,
        final_train_labels,
    )
    sl_fold_results_df.to_csv("results/sl_holdout_results.csv", index=False)

    # NG
    data = pd.read_csv(
        "/biogrid/BIOGRID-ORGANISM-Homo_sapiens-4.4.240.tab3.txt", sep="\t"
    )
    data = data[
        (data["Organism ID Interactor A"] == 9606)
        & (data["Organism ID Interactor B"] == 9606)
    ]
    data = data[data["Experimental System Type"] == "genetic"]

    selected_data = data[data["Experimental System"] == "Negative Genetic"]
    positive_pairs = []
    for i in range(selected_data.shape[0]):
        a = str(selected_data.iloc[i]["Entrez Gene Interactor A"])
        b = str(selected_data.iloc[i]["Entrez Gene Interactor B"])
        positive_pairs.append((min(a, b), max(a, b)))
    positive_pairs = set(positive_pairs)

    positive_pairs = {
        pair
        for pair in positive_pairs
        if pair[0] in reference_node2index and pair[1] in reference_node2index
    }
    print(len(positive_pairs))
    used_nodes = list(
        set([x for x, y in positive_pairs]).union(set([y for x, y in positive_pairs]))
    )

    holdout_nodes, splits, cv_nodes = helper.fold_split(used_nodes)
    (
        fold_splits,
        holdout_pairs,
        holdout_labels,
        final_train_pairs,
        final_train_labels,
    ) = helper.setup_data(positive_pairs, splits, holdout_nodes, cv_nodes)

    # save holdout and fold splits
    file_names = [
        "bin/ng_fold_splits_pairs.pkl",
        "bin/ng_fold_nodes.pkl",
        "bin/ng_holdout_pairs.pkl",
        "bin/ng_holdout_labels.pkl",
        "bin/ng_holdout_nodes.pkl",
    ]
    data_dicts = [fold_splits, splits, holdout_pairs, holdout_labels, holdout_nodes]

    for file_name, data_dict in zip(file_names, data_dicts):
        with open(file_name, "wb") as f:
            pickle.dump(data_dict, f)

    ng_fold_results_df, ng_holdout_results_df = helper.run_SVM(
        embeddings,
        reference_node2index,
        splits,
        fold_splits,
        holdout_pairs,
        holdout_labels,
        final_train_pairs,
        final_train_labels,
    )
    ng_holdout_results_df.to_csv(
        "results/ng_holdout_results.csv", index=False
    )

    # TF
    data = pd.read_csv("/data/tf_target.txt", sep="\t")
    tf_target_counts = data.groupby("TF").size()
    filtered_tfs = tf_target_counts[
        (tf_target_counts > 500) & (tf_target_counts < 1000)
    ].index
    print(len(filtered_tfs))
    selected_data = data[data["TF"].isin(filtered_tfs)]
    positive_pairs = []
    for i in range(selected_data.shape[0]):
        a = str(selected_data.iloc[i]["TF"])
        b = str(selected_data.iloc[i]["Target"])
        positive_pairs.append((min(a, b), max(a, b)))
    positive_pairs = set(positive_pairs)

    positive_pairs = {
        pair
        for pair in positive_pairs
        if pair[0] in reference_node2index and pair[1] in reference_node2index
    }
    print(len(positive_pairs))
    used_nodes = list(
        set([x for x, y in positive_pairs]).union(set([y for x, y in positive_pairs]))
    )

    holdout_nodes, splits, cv_nodes = helper.fold_split(used_nodes)
    (
        fold_splits,
        holdout_pairs,
        holdout_labels,
        final_train_pairs,
        final_train_labels,
    ) = helper.setup_data(positive_pairs, splits, holdout_nodes, cv_nodes)

    # save holdout and fold splits
    file_names = [
        "bin/tf_fold_splits_pairs.pkl",
        "bin/tf_fold_nodes.pkl",
        "bin/tf_holdout_pairs.pkl",
        "bin/tf_holdout_labels.pkl",
        "bin/tf_holdout_nodes.pkl",
    ]
    data_dicts = [fold_splits, splits, holdout_pairs, holdout_labels, holdout_nodes]

    for file_name, data_dict in zip(file_names, data_dicts):
        with open(file_name, "wb") as f:
            pickle.dump(data_dict, f)

    tf_fold_results_df, tf_holdout_results_df = helper.run_SVM(
        embeddings,
        reference_node2index,
        splits,
        fold_splits,
        holdout_pairs,
        holdout_labels,
        final_train_pairs,
        final_train_labels,
    )
    tf_holdout_results_df.to_csv(
        "results/tf_holdout_results.csv", index=False
    )

    # Compile data for plot
    meta_df = pd.read_csv("data/embed_meta.csv", index_col=0, encoding="utf-8")
    meta_df.index = meta_df.index.str.replace(r"\s+", "", regex=True)

    holdout_df = pd.concat(
        [sl_holdout_results_df, ng_holdout_results_df, tf_holdout_results_df],
        ignore_index=True,
    )
    holdout_auc = holdout_df.pivot(index="subfolder", columns="benchmark", values="AUC")
    holdout_auc = holdout_auc[["ssl", "ng", "tf"]]
    holdout_auc["average_auc"] = holdout_auc.mean(axis=1)
    holdout_auc = holdout_auc.sort_values(by="average_auc", ascending=False)
    holdout_auc = holdout_auc.drop(columns=["average_auc"])
    holdout_auc["data"] = meta_df.loc[holdout_auc.index, "Category"].values
    holdout_auc["algorithm"] = meta_df.loc[holdout_auc.index, "Method"].values
    holdout_auc["Dimensions"] = meta_df.loc[holdout_auc.index, "Dimensions"].values

    holdout_long = holdout_auc.reset_index().melt(
        id_vars=["subfolder", "data", "algorithm", "Dimensions"],
        value_vars=["sl", "ng", "tf"],
        var_name="Benchmark",
        value_name="AUC",
    )

    benchmark_colors = {"sl": "maroon", "ng": "navy", "tf": "#173317"}
    benchmark_markers = {"sl": "o", "ng": ">", "tf": "s"}

    plt.figure(figsize=(6.5, 8))
    for benchmark in holdout_long["Benchmark"].unique():
        benchmark_data = holdout_long[holdout_long["Benchmark"] == benchmark]
        sns.stripplot(
            data=benchmark_data,
            x="AUC",
            y="subfolder",
            color=benchmark_colors[benchmark],
            marker=benchmark_markers[benchmark],
            size=8,
            alpha=0.8,
            label=benchmark,
        )

    plt.title("", fontsize=16)
    plt.xlabel("AUC", fontsize=12)
    plt.ylabel("", fontsize=12)

    legend_elements = [
        Line2D([0], [0], color="maroon", marker="o", linestyle="", label="sl"),
        Line2D([0], [0], color="navy", marker=">", linestyle="", label="ng"),
        Line2D([0], [0], color="#173317", marker="s", linestyle="", label="tf"),
    ]

    plt.legend(handles=legend_elements, title="", loc="lower right", frameon=True)

    plt.tight_layout()
    plt.savefig(
        "results/plots/holdout_genepair_dot_plot.pdf", format="pdf", bbox_inches="tight"
    )

    plt.show()

    holdout_auc["Dimension_2"] = holdout_auc["Dimensions"].apply(lambda x: x[1])

    def calculate_anova_with_ratios(dependent_variable, data):
        print(f"\nANOVA for {dependent_variable}")
        model = smf.ols(
            f"{dependent_variable} ~ data + algorithm + Dimension_2", data=data
        ).fit()
        anova_table = sm.stats.anova_lm(model, typ=2)

        total_sum_sq = anova_table["sum_sq"].sum()
        anova_table["ratio"] = anova_table["sum_sq"] / total_sum_sq

        print(anova_table)
        return anova_table

    ssl_anova = calculate_anova_with_ratios("sl", holdout_auc)
    ng_anova = calculate_anova_with_ratios("ng", holdout_auc)
    tf_anova = calculate_anova_with_ratios("tf", holdout_auc)
