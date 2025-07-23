#!/usr/bin/env python3
import argparse
from pathlib import Path
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf

def calc_anova(dv, rhs, df, typ):
    model = smf.ols(f"{dv} ~ {rhs}", data=df).fit()
    tbl = sm.stats.anova_lm(model, typ=typ)
    total = tbl["sum_sq"].sum()
    tbl["ratio"] = tbl["sum_sq"] / total
    tbl.insert(0, "dependent_variable", dv)
    return tbl.reset_index().rename(columns={"index": "term"})

def main():
    ap = argparse.ArgumentParser(description="ANOVA with sum_sq ratios.")
    ap.add_argument("-i", "--input", required=True, help="Input CSV/TSV file of results.")
    ap.add_argument("--sep", help="Separator (default auto: ',' for .csv, '\\t' for .tsv).")
    ap.add_argument("-d", "--dv", required=True, nargs="+", help="Dependent variable column(s).")
    ap.add_argument("-t", "--terms", default="Method + Category + Dimension",
                    help='RHS of formula, e.g. "Method + Category + Dimension".')
    ap.add_argument("--typ", type=int, default=2, choices=[1,2,3],
                    help="Type of SS for ANOVA (1/2/3). Default 2.")
    ap.add_argument("-o", "--output", help="Output CSV path. Prints to stdout if omitted.")
    args = ap.parse_args()

    in_path = Path(args.input)
    sep = args.sep or ("\t" if in_path.suffix.lower() in {".tsv", ".tab"} else ",")
    df = pd.read_csv(in_path, sep=sep)

    out_frames = []
    for dv in args.dv:
        print(f"\n=== ANOVA for {dv} ===")
        res = calc_anova(dv, args.terms, df, args.typ)
        print(res)
        out_frames.append(res)

    out_df = pd.concat(out_frames, ignore_index=True)
    if args.output:
        out_df.to_csv(args.output, index=False)
    else:
        print("\n=== Combined results ===")
        print(out_df)

if __name__ == "__main__":
    # ./anova.py -i avg_results.tsv -d Y1 Y2 -t "Method + Category + Dimension" --typ 2 -o anova_out.csv
    main()
