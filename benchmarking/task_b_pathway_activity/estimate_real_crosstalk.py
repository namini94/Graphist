"""Closes the loop on the dose-response story (sim15) by asking the question directly on
REAL data, across MULTIPLE tissues: is a real phenotype-differentiated tissue
micro-environment actually associated with MORE multi-pathway crosstalk than a generic,
biologically-arbitrary region of the same tissue?

Generalized to run against any 10x Visium dataset with a real region/phenotype
annotation -- run three times this session:
  - lymph node: germinal center (GC) vs. non-GC (manual annotation from
    osmanbeyoglulab/STAN's own resources) -- GC reactions are a textbook coordinated
    multi-pathway immune process (chemokine + BCR + TCR + NF-kB signaling), literally
    the same program curated in simulate_phenotype_program.py's PHENOTYPE_PROGRAM.
  - BRCA-COMMOT (breast cancer Visium): Tumor vs. Healthy (iIMPACT manual annotation) --
    tumor regions are a textbook coordinated multi-pathway cancer process (EMT,
    hypoxia, proliferation, immune evasion all activate together).
  - DLPFC/Maynard (human cortex Visium, sample 151673): white matter (WM) vs. gray-
    matter layers (spatialLIBD expert annotation) -- myelination is a well-documented
    coordinated multi-gene program (oligodendrocyte differentiation + lipid/myelin
    biosynthesis).

Method (identical across all three, only the data/annotation loading differs):
1. Score all 30 panel pathways per spot via decoupleR GSVA (rank-based -- doesn't
   assume any pathway-pathway structure, so it can't manufacture crosstalk that isn't
   there) on real, library-size-normalized expression.
2. Crosstalk statistic for a given spot subset: mean |Pearson correlation| and count of
   BH-significant pairs, restricted to pathway pairs with DISJOINT gene membership
   (excludes pairs sharing genes, which would trivially correlate regardless of any real
   crosstalk -- same logic as simulate_phenotype_program.coordination_score, now applied
   to real inferred activity instead of synthetic ground truth).
3. Compare the real phenotype-group's statistic against a null built from many random
   same-sized subsets drawn from the comparison pool (typically the rest of the tissue,
   or a specific "generic" annotation value).

Usage:
  python estimate_real_crosstalk.py --h5 <filtered_feature_bc_matrix.h5> \
      --annot-file <annotation.csv> --annot-column <col> [--annot-sep ","] \
      --phenotype-values GC --background-values rest \
      --gmt <pathways.gmt> --out-dir <out> [--n-null 200] [--min-genes 5]
"""
import argparse
import itertools
import os

import decoupler as dc
import numpy as np
import pandas as pd
import scanpy as sc


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
    order = np.argsort(pvals)
    m = len(pvals)
    thresh = (np.arange(1, m + 1) / m) * fdr_alpha
    below = pvals[order] <= thresh
    n_sig = int(np.max(np.where(below)[0]) + 1) if below.any() else 0
    return {"mean_abs_r": float(rs.mean()), "n_sig_pairs": n_sig, "n_pairs": m}


def _normalize(counts: pd.DataFrame, panel_genes: list, source: str) -> pd.DataFrame:
    libsize = counts.sum(axis=1)
    keep = libsize > 0
    counts, libsize = counts[keep], libsize[keep]
    norm = np.log1p(counts.div(libsize, axis=0) * 1e4)
    print(f"Loaded {source}: {norm.shape[0]} spots x {norm.shape[1]} / {len(panel_genes)} panel genes present")
    return norm


def load_expression_h5(h5_path: str, panel_genes: list):
    """Loads a 10x filtered_feature_bc_matrix.h5, restricts to the panel's genes (as
    many as are present), library-size-normalizes (CP10K) and log1p's -- standard
    practice before GSVA. Returns (normalized expression, full original spot order --
    needed for positional annotation alignment, before any spots get dropped)."""
    adata = sc.read_10x_h5(h5_path)
    adata.var_names_make_unique()
    full_order = adata.obs_names.tolist()
    present = [g for g in panel_genes if g in adata.var_names]
    sub = adata[:, present]
    import scipy.sparse as sp
    X = sub.X.toarray() if sp.issparse(sub.X) else np.asarray(sub.X)
    counts = pd.DataFrame(X, index=sub.obs_names, columns=present)
    return _normalize(counts, panel_genes, h5_path), full_order


def load_expression_csv(counts_csv: str, panel_genes: list):
    """Loads a genes-x-spots counts.csv (this repo's convention, e.g.
    benchmarking/data/task_b/lymph_node_phenotype/counts.csv), restricts to the panel's
    genes present, transposes to spots x genes, normalizes as above. Returns (normalized
    expression, full original spot order)."""
    counts = pd.read_csv(counts_csv, index_col=0).T  # -> spots x genes
    full_order = counts.index.tolist()
    present = [g for g in panel_genes if g in counts.columns]
    return _normalize(counts[present], panel_genes, counts_csv), full_order


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5", default=None, help="10x filtered_feature_bc_matrix.h5")
    parser.add_argument("--counts-csv", default=None, help="genes x spots counts.csv (alternative to --h5)")
    parser.add_argument("--annot-file", required=True)
    parser.add_argument("--annot-sep", default=",")
    parser.add_argument("--annot-column", required=True)
    parser.add_argument("--annot-positional", action="store_true",
                         help="annotation file is indexed by row position (0-based order matching the "
                              "10x barcodes.tsv.gz / h5 spot order), not by spot barcode -- e.g. spatialLIBD's "
                              "DLPFC GT-labels.csv. Re-indexes the annotation onto the expression matrix's "
                              "own spot order positionally instead of matching on barcode strings.")
    parser.add_argument("--phenotype-values", required=True,
                         help="comma-separated annotation value(s) defining the real phenotype group")
    parser.add_argument("--background-values", default="rest",
                         help="comma-separated annotation value(s) defining the comparison pool, "
                              "or 'rest' for everything not in --phenotype-values")
    parser.add_argument("--gmt", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--n-null", type=int, default=200)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()
    rng = np.random.default_rng(args.seed)

    if bool(args.h5) == bool(args.counts_csv):
        raise ValueError("pass exactly one of --h5 or --counts-csv")
    gmt = read_gmt(args.gmt)
    panel_genes = sorted(set(g for gl in gmt.values() for g in gl))
    norm, full_spot_order = (load_expression_h5(args.h5, panel_genes) if args.h5
                              else load_expression_csv(args.counts_csv, panel_genes))

    net = dc.pp.read_gmt(args.gmt)
    scores, _ = dc.mt.gsva(data=norm, net=net, verbose=False)

    annot = pd.read_csv(args.annot_file, sep=args.annot_sep, index_col=0)
    if args.annot_positional:
        # Positional: row i of the annotation file corresponds to the i-th spot in the
        # ORIGINAL (pre-filtering) 10x spot order -- re-key onto those barcodes so
        # downstream isin()/index lookups work identically to the barcode-matched path.
        if len(annot) != len(full_spot_order):
            raise ValueError(f"positional annotation has {len(annot)} rows but the source has "
                              f"{len(full_spot_order)} spots -- can't align positionally")
        annot = annot.set_axis(full_spot_order)
    annot_col = annot[args.annot_column]
    pheno_values = set(args.phenotype_values.split(","))
    pheno_spots = [s for s in annot_col.index[annot_col.isin(pheno_values)] if s in scores.index]
    if args.background_values == "rest":
        bg_spots = [s for s in scores.index if s not in set(pheno_spots)]
    else:
        bg_values = set(args.background_values.split(","))
        bg_spots = [s for s in annot_col.index[annot_col.isin(bg_values)] if s in scores.index and s not in set(pheno_spots)]
    print(f"Phenotype spots ({sorted(pheno_values)}): {len(pheno_spots)}, "
          f"background pool ({args.background_values}): {len(bg_spots)}")

    pairs = disjoint_pairs(gmt)
    print(f"Disjoint-gene-membership pathway pairs available: {len(pairs)} / "
          f"{len(list(itertools.combinations(gmt.keys(), 2)))}")

    real_stat = crosstalk_stat(scores, pheno_spots, pairs)
    print(f"REAL phenotype spots -- mean|r|={real_stat['mean_abs_r']:.4f}  "
          f"significant pairs={real_stat['n_sig_pairs']}/{real_stat['n_pairs']}")

    null_mean_abs_r, null_n_sig = [], []
    for _ in range(args.n_null):
        # Bootstrap (with replacement): sound regardless of background-pool size. Sampling
        # WITHOUT replacement breaks down when the background pool is only modestly larger
        # than the phenotype group (e.g. Healthy n=485 vs. Tumor n=463 in BRCA-COMMOT) --
        # draws then overlap ~95%, collapsing the null's variance and producing an
        # unrealistically extreme z-score that isn't a fair test.
        draw = rng.choice(bg_spots, size=len(pheno_spots), replace=True).tolist()
        s = crosstalk_stat(scores, draw, pairs)
        null_mean_abs_r.append(s["mean_abs_r"])
        null_n_sig.append(s["n_sig_pairs"])
    null_mean_abs_r = np.array(null_mean_abs_r)
    null_n_sig = np.array(null_n_sig)

    z_mean_r = (real_stat["mean_abs_r"] - null_mean_abs_r.mean()) / (null_mean_abs_r.std() + 1e-12)
    pctile_mean_r = float((null_mean_abs_r < real_stat["mean_abs_r"]).mean() * 100)
    z_n_sig = (real_stat["n_sig_pairs"] - null_n_sig.mean()) / (null_n_sig.std() + 1e-12)
    pctile_n_sig = float((null_n_sig < real_stat["n_sig_pairs"]).mean() * 100)

    print(f"\nNull (random same-sized background subsets, n={args.n_null}):")
    print(f"  mean|r|: {null_mean_abs_r.mean():.4f} +/- {null_mean_abs_r.std():.4f}  "
          f"-> real phenotype z={z_mean_r:.2f}, percentile={pctile_mean_r:.1f}")
    print(f"  n_sig_pairs: {null_n_sig.mean():.1f} +/- {null_n_sig.std():.1f}  "
          f"-> real phenotype z={z_n_sig:.2f}, percentile={pctile_n_sig:.1f}")

    os.makedirs(args.out_dir, exist_ok=True)
    pd.DataFrame({"mean_abs_r": null_mean_abs_r, "n_sig_pairs": null_n_sig}).to_csv(
        os.path.join(args.out_dir, "null_distribution.csv"), index=False)
    pd.Series(real_stat).to_csv(os.path.join(args.out_dir, "real_gc_stat.csv"))
    print(f"\nWritten to {args.out_dir}")


if __name__ == "__main__":
    main()
