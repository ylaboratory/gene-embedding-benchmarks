#!/usr/bin/env python
import os

os.environ["OMP_NUM_THREADS"] = "10"
os.environ["MKL_NUM_THREADS"] = "10"
os.environ["NUMEXPR_NUM_THREADS"] = "10"
os.environ["VECLIB_MAXIMUM_THREADS"] = "10"
os.environ["OPENBLAS_NUM_THREADS"] = "10"
os.environ["BLIS_NUM_THREADS"] = "10"
import glob
import pickle
import time
import argparse

import numpy as np
import pandas as pd
from sklearn.svm import SVC
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import roc_auc_score, average_precision_score, make_scorer

C_VALUES = [0.1, 1, 10, 100, 1000]
N_JOBS = 3


def precision_at_k(y_true, y_scores, k=10):
    order = np.argsort(y_scores)[::-1]
    return np.mean(np.array(y_true)[order[:k]])


def filter_pairs_labels(pairs, labels, ref_genes):
    fp, fl, idx = [], [], []
    for i, ((g1, g2), lbl) in enumerate(zip(pairs, labels)):
        if g1 in ref_genes and g2 in ref_genes:
            fp.append((g1, g2))
            fl.append(lbl)
            idx.append(i)
    return fp, fl, idx


def build_features_sum(pairs, emb, g2i):
    feats = []
    for g1, g2 in pairs:
        feats.append(emb[g2i[g1]] + emb[g2i[g2]])
    return np.vstack(feats)


def build_features_product(pairs, emb, g2i):
    feats = []
    for g1, g2 in pairs:
        feats.append(emb[g2i[g1]] * emb[g2i[g2]])
    return np.vstack(feats)


def build_features_concat(pairs, emb, g2i):
    idx1 = [g2i[g1] for g1, g2 in pairs]
    idx2 = [g2i[g2] for g1, g2 in pairs]
    return np.hstack([emb[idx1], emb[idx2]])


def main():
    parser = argparse.ArgumentParser(description="Nested CV SVM on gene embeddings")
    parser.add_argument(
        "--subfolder",
        required=True,
        help="Path to embedding subfolder (with embedding CSV and genelist txt)",
    )
    parser.add_argument(
        "-o",
        "--operation",
        choices=["sum", "product", "concat"],
        default="sum",
        help="How to combine pair embeddings",
    )
    parser.add_argument(
        "-d", "--out-root", required=True, help="Directory to save results"
    )
    parser.add_argument("-s", "--suffix", help="Suffix for output CSV")
    parser.add_argument(
        "--cv-pkl", required=True, help="Path to nested CV splits pickle file"
    )
    args = parser.parse_args()

    subfolder = args.subfolder.rstrip("/")
    name = os.path.basename(subfolder)
    OUT_ROOT = args.out_root
    os.makedirs(OUT_ROOT, exist_ok=True)
    suffix = args.suffix or args.operation

    with open(args.cv_pkl, "rb") as f:
        cv_data = pickle.load(f)
    pairs, labels, cv_splits = cv_data["pairs"], cv_data["labels"], cv_data["cv_splits"]

    csvs = glob.glob(os.path.join(subfolder, "*.csv"))
    if not csvs:
        raise FileNotFoundError(f"No CSV files in {subfolder}")
    emb_vals = pd.read_csv(csvs[0], header=None).values

    txts = glob.glob(os.path.join(subfolder, "*.txt"))
    if len(txts) != 1:
        raise FileNotFoundError(
            f"Expected exactly 1 .txt gene-list in {subfolder}, found {len(txts)}"
        )
    with open(txts[0]) as f:
        ref_genes = [line.strip() for line in f if line.strip()]

    if args.operation == "sum":
        build_features = build_features_sum
    elif args.operation == "product":
        build_features = build_features_product
    else:
        build_features = build_features_concat

    g2i = {g: i for i, g in enumerate(ref_genes)}
    fp, fl, idx_master = filter_pairs_labels(pairs, labels, ref_genes)
    if not fp:
        raise ValueError("No pairs survive filtering for this embedding!")
    X = build_features(fp, emb_vals, g2i)
    y = np.array(fl)

    orig_pos = sum(l == 1 for l in labels)
    filtered_pos = sum(l == 1 for l in fl)
    print(f"Processing '{name}' | Operation: {args.operation}")
    print(f"Original positives: {orig_pos}, Filtered positives: {filtered_pos}")

    pr10_scorer = make_scorer(precision_at_k, needs_threshold=True, k=10)
    results = []

    # Outer CV loop
    for fold, splits in cv_splits.items():
        train_m, outer_m = splits["train_idx"], splits["test_idx"]
        train_idx = [i for i, m in enumerate(idx_master) if m in train_m]
        outer_idx = [i for i, m in enumerate(idx_master) if m in outer_m]

        X_train, y_train = X[train_idx], y[train_idx]
        X_outer, y_outer = X[outer_idx], y[outer_idx]

        # Map inner splits
        inner_splits_mapped = []
        train_master = [m for m in idx_master if m in train_m]
        full2train = {m: i for i, m in enumerate(train_master)}
        for tr_m, val_m in splits["inner_splits"]:
            tr_idx = [full2train[m] for m in tr_m if m in full2train]
            val_idx = [full2train[m] for m in val_m if m in full2train]
            if tr_idx and val_idx:
                inner_splits_mapped.append((tr_idx, val_idx))

        print(f"Fold {fold}: train={len(train_idx)}, outer={len(outer_idx)}")

        # Inner CV
        t0_inner = time.time()
        grid = GridSearchCV(
            SVC(class_weight="balanced", probability=False),
            param_grid={"C": C_VALUES},
            scoring={
                "AUC": "roc_auc",
                "AUPRC": "average_precision",
                "PR@10": pr10_scorer,
            },
            refit="AUC",
            cv=inner_splits_mapped,
            n_jobs=N_JOBS,
            return_train_score=False,
        )
        grid.fit(X_train, y_train)
        inner_time = time.time() - t0_inner

        best_i = grid.best_index_
        cvres = grid.cv_results_
        inner_auc = grid.best_score_
        inner_auprc = cvres["mean_test_AUPRC"][best_i]
        inner_pr10 = cvres["mean_test_PR@10"][best_i]

        # Outer evaluation
        t0_outer = time.time()
        outer_scores = grid.decision_function(X_outer)
        outer_auc = roc_auc_score(y_outer, outer_scores)
        outer_auprc = average_precision_score(y_outer, outer_scores)
        outer_pr10 = precision_at_k(y_outer, outer_scores, k=10)
        outer_time = time.time() - t0_outer

        results.append(
            {
                "fold": fold,
                "best_C": grid.best_params_["C"],
                "inner_AUC": inner_auc,
                "inner_AUPRC": inner_auprc,
                "inner_PR@10": inner_pr10,
                "inner_time": inner_time,
                "outer_AUC": outer_auc,
                "outer_AUPRC": outer_auprc,
                "outer_PR@10": outer_pr10,
                "outer_time": outer_time,
            }
        )

    df = pd.DataFrame(results)
    df["orig_pos_pairs"] = orig_pos
    df["filtered_pos_pairs"] = filtered_pos

    summary = {
        "fold": "average",
        "best_C": df["best_C"].mode()[0],
        "inner_AUC": df["inner_AUC"].mean(),
        "inner_AUPRC": df["inner_AUPRC"].mean(),
        "inner_PR@10": df["inner_PR@10"].mean(),
        "inner_time": df["inner_time"].mean(),
        "outer_AUC": df["outer_AUC"].mean(),
        "outer_AUPRC": df["outer_AUPRC"].mean(),
        "outer_PR@10": df["outer_PR@10"].mean(),
        "outer_time": df["outer_time"].mean(),
        "orig_pos_pairs": df["orig_pos_pairs"].mean(),
        "filtered_pos_pairs": df["filtered_pos_pairs"].mean(),
    }
    df = pd.concat([df, pd.DataFrame([summary])], ignore_index=True)

    out_csv = os.path.join(OUT_ROOT, f"{name}_{suffix}.csv")
    df.to_csv(out_csv, index=False)

    print("\nResults:")
    print(df.to_string(index=False))
    print(f"\nSaved to {out_csv}")


if __name__ == "__main__":
    main()
