import pandas as pd
import mygene


def process_embedding_files(
    file_names, scopes="uniprot", base_path="data/embeddings/original"
):
    mg = mygene.MyGeneInfo()

    processed_dfs = {}
    entrez_sets = []

    for file in file_names:
        if file.endswith(".csv"):
            filename = base_path + file
            df = pd.read_csv(filename)
            print(f"Processing {file}...")

            uniprot_ids = df["Entry"].tolist()
            mg_results = mg.querymany(
                uniprot_ids, scopes=scopes, fields="entrezgene", species="human"
            )
            mg_df = pd.DataFrame(mg_results)

            df_with_entrez = df.merge(
                mg_df[["query", "entrezgene"]], left_on="Entry", right_on="query"
            )
            df_with_entrez.dropna(subset=["entrezgene"], inplace=True)
            df_with_entrez.drop(columns=["query", "Entry"], inplace=True)
            df_with_entrez.reset_index(inplace=True)

            df_grouped = (
                df_with_entrez.groupby("entrezgene", sort=False)
                .agg(
                    {
                        "index": "min",
                        **{
                            col: "mean"
                            for col in df_with_entrez.columns
                            if col not in ["entrezgene", "index"]
                        },
                    }
                )
                .reset_index()
            )
            df_grouped.sort_values("index", inplace=True)
            df_grouped.drop(columns=["index"], inplace=True)
            print(df_grouped.shape)

            processed_dfs[file] = df_grouped
            entrez_set = set(df_grouped["entrezgene"])
            print(f"Initial size: {len(entrez_set)}")
            entrez_sets.append(entrez_set)
        else:
            print(f"{file} is not a CSV...")

    print("Processing complete.")
    return processed_dfs, entrez_sets
