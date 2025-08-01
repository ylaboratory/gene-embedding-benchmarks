import os
os.environ["OMP_NUM_THREADS"] = "10"  
os.environ["MKL_NUM_THREADS"] = "10" 
os.environ["NUMEXPR_NUM_THREADS"] = "10"  
os.environ["VECLIB_MAXIMUM_THREADS"] = "10"  
os.environ["OPENBLAS_NUM_THREADS"] = "10" 
import pandas as pd
import numpy as np
import os
import glob
from sklearn.metrics import roc_auc_score, average_precision_score, precision_score
from sklearn.svm import SVC
import time
import pickle
import argparse


def precision_at_k(y_true, y_scores, k=10):
    if len(y_scores) > k:
        idx = np.argsort(y_scores)[-k:][::-1]
    else:
        idx = np.argsort(y_scores)[::-1]
    return np.mean(np.array(y_true)[idx])


def get_xy(df, emb_df):
    genes = df['gene'].values
    present = emb_df.index.intersection(genes)
    filtered = df[df['gene'].isin(present)]
    X = emb_df.loc[filtered['gene'].values]
    y = filtered['result'].values
    return X, y

C_values = [0.1, 1.0, 10.0, 100.0, 1000.0]


def main():
    parser = argparse.ArgumentParser(description="task-term SVM evaluation with nested CV")
    parser.add_argument('--subfolder', required=True,
                        help="Path to embedding subfolder containing CSV and gene-list .txt")
    parser.add_argument('--cv-fold1-pkl', required=True,
                        help="Path to task_cv_fold1_dict_all.pkl")
    parser.add_argument('--cv-fold2-pkl', required=True,
                        help="Path to task_cv_fold2_dict_all.pkl")
    parser.add_argument('--cv-fold3-pkl', required=True,
                        help="Path to task_cv_fold3_dict_all.pkl")
    parser.add_argument('--holdout-pkl', required=True,
                        help="Path to task_holdout_dict.pkl")
    parser.add_argument('-d', '--out-root', required=True,
                        help="Output directory to save results")
    args = parser.parse_args()

    subfolder = args.subfolder.rstrip('/')
    name = os.path.basename(subfolder)
    os.makedirs(args.out_root, exist_ok=True)

    csvs = glob.glob(os.path.join(subfolder, '*.csv'))
    if not csvs:
        raise FileNotFoundError(f"No CSV in {subfolder}")
    emb = pd.read_csv(csvs[0], header=None).values
    txts = glob.glob(os.path.join(subfolder, '*.txt'))
    if len(txts) != 1:
        raise FileNotFoundError(f"Expected one .txt gene list in {subfolder}, found {len(txts)}")
    with open(txts[0]) as f:
        genes = [l.strip() for l in f if l.strip()]
    embedding_df = pd.DataFrame(emb, index=genes)

    with open(args.cv_fold1_pkl, 'rb') as f:
        cv_fold1_dict = pickle.load(f)
    with open(args.cv_fold2_pkl, 'rb') as f:
        cv_fold2_dict = pickle.load(f)
    with open(args.cv_fold3_pkl, 'rb') as f:
        cv_fold3_dict = pickle.load(f)
    with open(args.holdout_pkl, 'rb') as f:
        holdout_dict = pickle.load(f)

    all_fold_results = {}
    all_holdout_results = {}

    for task_term in holdout_dict.keys():
        print(task_term, end =",", flush=True)
        
        fold1_df = cv_fold1_dict[task_term]
        fold2_df = cv_fold2_dict[task_term]
        fold3_df = cv_fold3_dict[task_term]
        holdout_df = holdout_dict[task_term]
        
        fold_results_records = []
        holdout_results_records = []
        
        
        fold_splits = [
            (fold1_df, pd.concat([fold2_df, fold3_df], ignore_index=True)),
            (fold2_df, pd.concat([fold1_df, fold3_df], ignore_index=True)),
            (fold3_df, pd.concat([fold1_df, fold2_df], ignore_index=True))
        ]
        
        c_performance_records = []
        
        for C_val in C_values:
            C = C_val
            fold_results = []
            fold_times = []
            
            for i, (test_df, train_df) in enumerate(fold_splits, start=1):
                # Start timer for this fold
                fold_start_time = time.time()
                
                X_train, y_train = get_xy(train_df, embedding_df)
                X_test, y_test = get_xy(test_df, embedding_df)

                if len(X_train) == 0 or len(X_test) == 0:
                    fold_results.append((np.nan, np.nan))
                    fold_times.append(np.nan)
                else:
                    clf = SVC(C=C_val, class_weight='balanced', probability=True)
                    if len(np.unique(y_train)) < 2:
                        print(f"Skipping {task_term} fold {i}: only one class present.")
                        continue

                    clf.fit(X_train, y_train)

                    y_prob = clf.decision_function(X_test)

                    auc = roc_auc_score(y_test, y_prob) if len(np.unique(y_test)) > 1 else np.nan
                    auprc = average_precision_score(y_test, y_prob)
                    pr10 = precision_at_k(y_test, y_prob, k=10)
                    
                    # End timer for this fold
                    fold_end_time = time.time()
                    fold_duration = fold_end_time - fold_start_time
                    fold_times.append(fold_duration)



                    fold_results.append((auc, auprc, pr10))
            
            avg_auc = np.mean([r[0] for r in fold_results])
            avg_auprc = np.mean([r[1] for r in fold_results])
            avg_pr10 = np.mean([r[2] for r in fold_results])
            avg_fold_time = np.mean(fold_times)  
            
            c_performance_records.append({
                'C': C_val,
                'AUC': avg_auc,
                'AUPRC': avg_auprc,
                'PR@10': avg_pr10,
                'avg_fold_time': avg_fold_time  
            })

        
        best_C_entry = max(c_performance_records, key=lambda x: x['AUC'])
        best_C = best_C_entry['C']
        best_avg_auc = best_C_entry['AUC']
        best_avg_auprc = best_C_entry['AUPRC']
        best_avg_pr10 = best_C_entry['PR@10']
        best_avg_fold_time = best_C_entry['avg_fold_time']


        fold_results_records.append({
            'subfolder': subfolder,
            'C': best_C,
            'AUC': best_avg_auc,
            'AUPRC': best_avg_auprc,
            'PR@10': best_avg_pr10,
            'avg_fold_time': best_avg_fold_time
        })
        
        df_all__ = pd.concat([fold1_df, fold2_df, fold3_df], ignore_index=True)
        unique_before = df_all__['gene'].nunique()
        print(f"Unique genes before filter: {unique_before}")
        genes_in_emb__ = embedding_df.index.intersection(df_all__['gene'])
        df_filtered__ = df_all__[df_all__['gene'].isin(genes_in_emb__)]
        unique_after = df_filtered__['gene'].nunique()
        print(f"Unique genes in embed: {unique_after}")

        X_train_all, y_train_all = get_xy(pd.concat([fold1_df, fold2_df, fold3_df], ignore_index=True), embedding_df)
        X_holdout, y_holdout = get_xy(holdout_df, embedding_df)
        
        holdout_start_time = time.time()
        
        if len(X_train_all) == 0 or len(X_holdout) == 0:
            holdout_auc = np.nan
            holdout_auprc = np.nan
        else:
            clf = SVC(C=best_C, class_weight='balanced', probability=True)
            if len(np.unique(y_train_all)) < 2:
                print(f"Skipping {task_term}: only one class present.")
                continue

            clf.fit(X_train_all, y_train_all)
            y_prob_holdout = clf.decision_function(X_holdout)

            holdout_auc = (roc_auc_score(y_holdout, y_prob_holdout) 
                        if len(np.unique(y_holdout)) > 1 else np.nan)
            holdout_auprc = average_precision_score(y_holdout, y_prob_holdout)
            holdout_pr10 = precision_at_k(y_holdout, y_prob_holdout, k=10)

        
        holdout_end_time = time.time()
        holdout_duration = holdout_end_time - holdout_start_time
        
        
        print(f"Holdout performance for {subfolder} with best C={best_C}:")
        print(f"Holdout AUC: {holdout_auc:.4f}")
        print(f"Holdout AUPRC: {holdout_auprc:.4f}")
        print(f"Holdout PR@10: {holdout_pr10:.4f}")

        print(f"Holdout Duration: {holdout_duration:.2f}s", flush = True)

        holdout_results_records.append({
            'subfolder': subfolder,
            'C': best_C,
            'AUC': holdout_auc,
            'AUPRC': holdout_auprc,
            'PR@10': holdout_pr10,
            'holdout_time': holdout_duration,
            'training gene in emb': unique_after,
            'traning gene all': unique_before
        })

        
        fold_results_df = pd.DataFrame(fold_results_records).set_index('subfolder')
        holdout_results_df = pd.DataFrame(holdout_results_records).set_index('subfolder')
        
        all_fold_results[task_term] = fold_results_df
        all_holdout_results[task_term] = holdout_results_df

    # save results
    base = os.path.basename(subfolder)
    out1 = os.path.join(args.out_root, f"{base}_fold_results.pkl")
    out2 = os.path.join(args.out_root, f"{base}_holdout_results.pkl")
    with open(out1, 'wb') as f:
        pickle.dump(all_fold_results, f)
    with open(out2, 'wb') as f:
        pickle.dump(all_holdout_results, f)
    print(f"Saved fold results to {out1}\nSaved holdout results to {out2}")

if __name__ == '__main__':
    main()
