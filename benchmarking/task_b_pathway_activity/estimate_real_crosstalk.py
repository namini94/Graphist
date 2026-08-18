"""Closes the loop on the dose-response story (sim15) by asking the question directly on
REAL data: is a real phenotype-differentiated tissue micro-environment actually
associated with MORE multi-pathway crosstalk than a generic, biologically-arbitrary
region of the same tissue?

Uses the 10x human lymph node Visium data (same dataset as sim6/16_realistic) with a
real, established phenotype axis: germinal center (GC) vs. non-GC spots (manual
annotation from `osmanbeyoglulab/STAN`'s own resources -- the same dataset/annotation
STAN itself validates against). GC reactions are a textbook coordinated multi-pathway
immune process -- chemokine-driven B-cell recruitment, BCR signaling, T-cell help, NF-kB
activation -- literally the same immune program curated in simulate_phenotype_program.py's
PHENOTYPE_PROGRAM, making this a direct, not just analogous, real-data test.

Method:
1. Score all 30 panel pathways per spot on REAL, normalized lymph node expression via
   decoupleR (GSVA -- rank-based, doesn't assume any pathway-pathway structure, so it
   can't manufacture crosstalk that isn't there).
2. Crosstalk statistic for a given spot subset: mean |Pearson correlation| across all
   pathway PAIRS WITH DISJOINT gene membership (excludes pairs sharing genes, which
   would trivially correlate regardless of any real crosstalk) -- same logic as
   simulate_phenotype_program.coordination_score, now applied to real inferred activity
   instead of synthetic ground truth.
3. Compare the REAL GC-spot subset's crosstalk statistic against a null distribution
   built from many random same-sized subsets of non-GC spots -- "how much more
   coordinated is the real phenotype micro-environment than an arbitrary same-sized
   piece of the same tissue?"

Usage: python estimate_real_crosstalk.py --data-dir <lymph_node_phenotype> \
    --gc-annot <manual_GC_annot.csv> --gmt <pathways.gmt> --out-dir <out> [--n-null 200]
"""
import argparse
import itertools
import os
import sys

import decoupler as dc
import numpy as np
import pandas as pd


def read_gmt(path):
    d = {}
    with open(path) as f:
        for line in f:
            parts = line.strip().split("\t")
            d[parts[0]] = parts[2:]
    return d


def disjoint_pairs(gmt: dict) -> list:
    pairs = []
    for p1, p2 in itertools.combinations(gmt.keys(), 2):
        if set(gmt[p1]).isdisjoint(gmt[p2]):
            pairs.append((p1, p2))
    return pairs


def crosstalk_stat(scores: pd.DataFrame, spot_ids: list, pairs: list, fdr_alpha: float = 0.05) -> dict:
    from scipy.stats import pearsonr
    sub = scores.loc[spot_ids]
    rs, pvals = [], []
    for p1, p2 in pairs:
        r, p = pearsonr(sub[p1].values, sub[p2].values)
        rs.append(abs(r))
        pvals.append(p)
    rs, pvals = np.array(rs), np.array(pvals)
    # BH correction (same convention as evaluate.py's benjamini_hochberg)
    order = np.argsort(pvals)
    m = len(pvals)
    thresh = (np.arange(1, m + 1) / m) * fdr_alpha
    below = pvals[order] <= thresh
    n_sig = int(np.max(np.where(below)[0]) + 1) if below.any() else 0
    return {"mean_abs_r": float(rs.mean()), "n_sig_pairs": n_sig, "n_pairs": m}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", required=True, help="dir with counts.csv (genes x spots, raw)")
    parser.add_argument("--gc-annot", required=True)
    parser.add_argument("--gmt", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--n-null", type=int, default=200)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()
    rng = np.random.default_rng(args.seed)

    counts = pd.read_csv(os.path.join(args.data_dir, "counts.csv"), index_col=0)  # genes x spots
    expr = counts.T  # spots x genes, decoupleR's expected orientation
    libsize = expr.sum(axis=1)
    norm = np.log1p(expr.div(libsize, axis=0) * 1e4)  # CP10K + log1p, standard practice

    gmt = read_gmt(args.gmt)
    net = dc.pp.read_gmt(args.gmt)
    scores, _ = dc.mt.gsva(data=norm, net=net, verbose=False)

    gc_annot = pd.read_csv(args.gc_annot, index_col=0)
    gc_spots = gc_annot.index[gc_annot.iloc[:, 0] == "GC"].tolist()
    gc_spots = [s for s in gc_spots if s in scores.index]
    non_gc_spots = [s for s in scores.index if s not in set(gc_spots)]
    print(f"GC spots: {len(gc_spots)}, non-GC spots: {len(non_gc_spots)}")

    pairs = disjoint_pairs(gmt)
    print(f"Disjoint-gene-membership pathway pairs available: {len(pairs)} / {len(list(itertools.combinations(gmt.keys(), 2)))}")

    real_stat = crosstalk_stat(scores, gc_spots, pairs)
    print(f"REAL GC spots  -- mean|r|={real_stat['mean_abs_r']:.4f}  "
          f"significant pairs={real_stat['n_sig_pairs']}/{real_stat['n_pairs']}")

    null_mean_abs_r, null_n_sig = [], []
    for i in range(args.n_null):
        draw = rng.choice(non_gc_spots, size=len(gc_spots), replace=False).tolist()
        s = crosstalk_stat(scores, draw, pairs)
        null_mean_abs_r.append(s["mean_abs_r"])
        null_n_sig.append(s["n_sig_pairs"])
    null_mean_abs_r = np.array(null_mean_abs_r)
    null_n_sig = np.array(null_n_sig)

    z_mean_r = (real_stat["mean_abs_r"] - null_mean_abs_r.mean()) / (null_mean_abs_r.std() + 1e-12)
    pctile_mean_r = float((null_mean_abs_r < real_stat["mean_abs_r"]).mean() * 100)
    z_n_sig = (real_stat["n_sig_pairs"] - null_n_sig.mean()) / (null_n_sig.std() + 1e-12)
    pctile_n_sig = float((null_n_sig < real_stat["n_sig_pairs"]).mean() * 100)

    print(f"\nNull (random same-sized non-GC subsets, n={args.n_null}):")
    print(f"  mean|r|: {null_mean_abs_r.mean():.4f} +/- {null_mean_abs_r.std():.4f}  "
          f"-> real GC z={z_mean_r:.2f}, percentile={pctile_mean_r:.1f}")
    print(f"  n_sig_pairs: {null_n_sig.mean():.1f} +/- {null_n_sig.std():.1f}  "
          f"-> real GC z={z_n_sig:.2f}, percentile={pctile_n_sig:.1f}")

    os.makedirs(args.out_dir, exist_ok=True)
    pd.DataFrame({"mean_abs_r": null_mean_abs_r, "n_sig_pairs": null_n_sig}).to_csv(
        os.path.join(args.out_dir, "null_distribution.csv"), index=False)
    pd.Series(real_stat).to_csv(os.path.join(args.out_dir, "real_gc_stat.csv"))
    scores.to_csv(os.path.join(args.out_dir, "pathway_scores.csv"))
    gc_annot.loc[scores.index].to_csv(os.path.join(args.out_dir, "gc_annotation_matched.csv"))
    print(f"\nWritten to {args.out_dir}")


if __name__ == "__main__":
    main()
