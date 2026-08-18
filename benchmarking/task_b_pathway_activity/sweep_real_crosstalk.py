"""Exhaustive version of estimate_real_crosstalk.py: instead of one hand-picked
phenotype-vs-background pair per dataset, tests EVERY pairwise combination of a
dataset's annotation categories -- pre-empting the reasonable reviewer objection that
GC/Tumor/gray-matter were cherry-picked to produce the desired result.

For each pair of categories, the smaller group is always the one tested directly
(the well-powered direction -- see estimate_real_crosstalk.py's docstring on why
bootstrapping a small pool up to a much larger target corrupts the significance-count
metric) against a bootstrap null built from the larger group.

GSVA is computed once per (dataset, panel) and reused across all pairs, rather than
re-invoking the whole pipeline per comparison.

Usage: python sweep_real_crosstalk.py --h5 <...> --annot-file <...> --annot-column <col> \
    [--annot-sep ","] [--annot-positional] --gmt <...> --out-csv <...> \
    [--min-group-size 50] [--n-null 100] [--seed 0]
"""
import argparse
import itertools
import os

import decoupler as dc
import numpy as np
import pandas as pd

from estimate_real_crosstalk import (
    read_gmt, disjoint_pairs, crosstalk_stat, load_expression_h5, load_expression_csv,
)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5", default=None)
    parser.add_argument("--counts-csv", default=None)
    parser.add_argument("--annot-file", required=True)
    parser.add_argument("--annot-sep", default=",")
    parser.add_argument("--annot-column", required=True)
    parser.add_argument("--annot-positional", action="store_true")
    parser.add_argument("--gmt", required=True)
    parser.add_argument("--out-csv", required=True)
    parser.add_argument("--min-group-size", type=int, default=50,
                         help="skip categories with fewer spots than this -- too small to give a "
                              "reliable per-pathway-pair correlation estimate")
    parser.add_argument("--n-null", type=int, default=100,
                         help="fewer than the headline analyses' 200, since this runs many "
                              "comparisons -- still enough for a robustness check, not the primary claim")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--dataset-label", default="")
    args = parser.parse_args()
    rng = np.random.default_rng(args.seed)

    gmt = read_gmt(args.gmt)
    panel_genes = sorted(set(g for gl in gmt.values() for g in gl))
    norm, full_spot_order = (load_expression_h5(args.h5, panel_genes) if args.h5
                              else load_expression_csv(args.counts_csv, panel_genes))
    net = dc.pp.read_gmt(args.gmt)
    scores, _ = dc.mt.gsva(data=norm, net=net, verbose=False)

    annot = pd.read_csv(args.annot_file, sep=args.annot_sep, index_col=0)
    if args.annot_positional:
        if len(annot) != len(full_spot_order):
            raise ValueError("positional annotation length mismatch")
        annot = annot.set_axis(full_spot_order)
    annot_col = annot[args.annot_column].reindex(scores.index)

    counts = annot_col.value_counts()
    categories = [c for c in counts.index if counts[c] >= args.min_group_size]
    print(f"Categories with >= {args.min_group_size} spots: "
          + ", ".join(f"{c} (n={counts[c]})" for c in categories))

    pairs = disjoint_pairs(gmt)
    print(f"Disjoint pathway pairs: {len(pairs)}")

    rows = []
    for cat_a, cat_b in itertools.combinations(categories, 2):
        spots_a = annot_col.index[annot_col == cat_a].tolist()
        spots_b = annot_col.index[annot_col == cat_b].tolist()
        # well-powered direction: test the smaller group against the larger as background
        if len(spots_a) <= len(spots_b):
            small, small_name, large, large_name = spots_a, cat_a, spots_b, cat_b
        else:
            small, small_name, large, large_name = spots_b, cat_b, spots_a, cat_a

        real_stat = crosstalk_stat(scores, small, pairs)
        null_r, null_n = [], []
        for _ in range(args.n_null):
            draw = rng.choice(large, size=len(small), replace=True).tolist()
            s = crosstalk_stat(scores, draw, pairs)
            null_r.append(s["mean_abs_r"])
            null_n.append(s["n_sig_pairs"])
        null_r, null_n = np.array(null_r), np.array(null_n)
        z_r = (real_stat["mean_abs_r"] - null_r.mean()) / (null_r.std() + 1e-12)
        z_n = (real_stat["n_sig_pairs"] - null_n.mean()) / (null_n.std() + 1e-12)

        row = {
            "dataset": args.dataset_label, "tested_group": small_name, "background_group": large_name,
            "tested_n": len(small), "background_n": len(large),
            "mean_abs_r": real_stat["mean_abs_r"], "null_mean_abs_r": null_r.mean(), "null_sd_abs_r": null_r.std(),
            "z_strength": z_r,
            "n_sig_pairs": real_stat["n_sig_pairs"], "null_mean_n_sig": null_n.mean(), "null_sd_n_sig": null_n.std(),
            "z_breadth": z_n,
        }
        rows.append(row)
        print(f"  {small_name} (n={len(small)}) vs {large_name} (n={len(large)}): "
              f"z_strength={z_r:+.2f}  z_breadth={z_n:+.2f}")

    df = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(args.out_csv), exist_ok=True)
    df.to_csv(args.out_csv, index=False)
    print(f"\nWritten {len(df)} pairwise comparisons to {args.out_csv}")


if __name__ == "__main__":
    main()
