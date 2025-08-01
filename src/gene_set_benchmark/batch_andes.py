import os
import subprocess
from datetime import datetime
import argparse
from pathlib import Path

# change file name/path (that are commented out) based on disease-tissue vs kegg-go task

parser = argparse.ArgumentParser()

parser.add_argument("--embed-path", default="data/embeddings/intersect", type=Path)
parser.add_argument("--andes-script", default="ANDES/src/andes.py", type=Path)
parser.add_argument("--geneset1", default="data/gmt/KEGG_CPDB.gmt", type=Path)
parser.add_argument(
    "--geneset2", default="data/gmt/hsa_low_eval_BP_propagated.gmt", type=Path
)
parser.add_argument("--out-dir", default="results/andes_out", type=Path)
parser.add_argument("--status-log", default="andes_status.txt", type=Path)
parser.add_argument(
    "--matched-pairs", default="data/matched_pairs/matched_pair_kegg_go.txt", type=Path
)

args = parser.parse_args()

embed_path = args.embed_path
andes_script = args.andes_script
geneset1 = args.geneset1
geneset2 = args.geneset2
andes_out_path = args.out_dir
status_log = args.status_log
matched_pairs_file = args.matched_pairs

# embed_path = "data/embeddings/intersect"
# andes_script = "ANDES/src/andes.py"
# geneset1 = "data/gmt/KEGG_CPDB.gmt"  # 'data/gmt/bto_specific.gmt'
# geneset2 = "data/gmt/hsa_low_eval_BP_propagated.gmt"  # 'data/gmt/omim_entrez.gmt'
# andes_out_path = "results/andes_out"
# status_log = "andes_status.txt"
# matched_pairs_file = "data/matched_pairs/matched_pair_kegg_go.txt"  # data/matched_pairs/matched_list_bto_doid.txt

if not os.path.exists(andes_out_path):
    os.makedirs(andes_out_path)

with open(status_log, "a") as log_file:
    log_file.write(f"{datetime.now()}: new run of ANDES tests\n")

    subdirs = [
        d for d in os.listdir(embed_path) if os.path.isdir(os.path.join(embed_path, d))
    ]
    print(subdirs)

    for subdir in subdirs:
        print(subdir)
        subdir_path = os.path.join(embed_path, subdir)
        files = os.listdir(subdir_path)
        csv_files = [f for f in files if f.endswith(".csv")]
        txt_files = [f for f in files if f.endswith(".txt")]

        if len(csv_files) == 0:
            log_file.write(f"{datetime.now()}: No CSV file found in {subdir_path}\n")
            log_file.flush()
            continue
        if len(txt_files) == 0:
            log_file.write(f"{datetime.now()}: No TXT file found in {subdir_path}\n")
            log_file.flush()
            continue

        emb_file = os.path.join(subdir_path, csv_files[0])
        genelist_file = os.path.join(subdir_path, txt_files[0])

        output_filename = f"{subdir}_andesout.csv"
        output_file = os.path.join(andes_out_path, output_filename)

        start_time = datetime.now()
        log_file.write(f"{start_time}: Starting ANDES on {subdir}\n")

        command = [
            "python",
            andes_script,
            "--emb",
            emb_file,
            "--genelist",
            genelist_file,
            "--geneset1",
            geneset1,
            "--geneset2",
            geneset2,
            "--out",
            output_file,
            "-n",
            "25",
            "--matchedpair",
            matched_pairs_file,
        ]

        log_file.write(f"Using this command: {' '.join(command)}\n")
        log_file.flush()

        try:
            subprocess.run(command, check=True)
            end_time = datetime.now()
            log_file.write(f"{end_time}: Finished ANDES on {subdir}\n")
            log_file.flush()
        except subprocess.CalledProcessError as e:
            end_time = datetime.now()
            log_file.write(f"{end_time}: Error running ANDES on {subdir}: {e}\n")
            log_file.flush()

    final_time = datetime.now()
    log_file.write(f"{final_time}: All processing completed.\n")
    log_file.flush()

print("done")
