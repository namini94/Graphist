"""Task B "realistic-noise" scenario: known pathway activity injected on top of a
real-data-grounded count-noise backbone, instead of the additive-Gaussian-noise model
used by sim1-sim5.

Motivation (user's own concern): the hand-rolled Gaussian-noise simulator in
simulate_pathway_activity.py has no guarantee its noise characteristics look anything like
real ST data. This scenario addresses that directly: the per-spot, per-gene (mean,
dispersion) backbone comes from fit_realistic_backbone.R, which fits scDesign3's
negative-binomial marginal model (spatially-smooth mean via a GP spline, gene-specific
dispersion) to REAL 10x human lymph node Visium data, then evaluates it at our synthetic
30x30 grid coordinates. This script:

1. Reuses that realistic (mu, sigma) backbone as-is -- no change to the real-data-fitted
   noise structure.
2. Re-simulates a known per-spot, per-pathway "true activity" ground truth (same
   spatially-smooth-field + group-demeaning + injected-DE-shift logic as
   simulate_pathway_activity.simulate_true_activity, for consistency with sim1-sim5),
   using the SAME 30x30 grid / A-B group split the backbone was evaluated at (so the known
   effect lines up spatially with the realistic mean surface).
3. Injects that activity multiplicatively on the log-mean (standard NB-GLM log link):
   mu_final = realistic_mu * exp(effect_scale * activity @ weights.T)
   then samples counts ~ NegativeBinomial(mean=mu_final, dispersion=realistic_sigma) --
   dispersion (gamlss NBI convention, Var = mu + sigma*mu^2) is converted to numpy's
   (n, p) parameterization via n = 1/sigma, p = n/(n+mu).
4. Writes log1p(counts) as st_expression.csv -- the same interface run_graphist.py /
   run_stan.py / run_decoupler.py already consume, so every existing baseline script works
   unchanged on this scenario.

The pathway panel is necessarily a SUBSET of sim1's 30 pathways: only pathway member genes
that were actually detected in the real lymph node data (600 of the panel's 732 genes) can
get a realistic backbone value. Pathways left with too few surviving genes are dropped.

Usage: python simulate_realistic_scenario.py --backbone-dir <sim6_realistic> \
    --template-dir <sim1> --out-dir <out>
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from simulate_pathway_activity import read_gmt, simulate_true_activity, write_gmt  # noqa: E402


def subset_pathway_mask(gmt_path: str, genes_available: list, min_genes: int):
    """Restrict a pathway .gmt to genes present in the realistic backbone; drop pathways
    left with fewer than min_genes surviving members. Returns (mask [genes x pathways]
    aligned to genes_available's order, pathway_names, subsetted_gmt_dict)."""
    all_pathways = read_gmt(gmt_path)
    available = set(genes_available)
    subsetted = {}
    for name, members in all_pathways.items():
        kept = [g for g in members if g in available]
        if len(kept) >= min_genes:
            subsetted[name] = kept

    pathway_names = list(subsetted.keys())
    gene_idx = {g: i for i, g in enumerate(genes_available)}
    mask = np.zeros((len(genes_available), len(pathway_names)), dtype=float)
    for j, name in enumerate(pathway_names):
        for g in subsetted[name]:
            mask[gene_idx[g], j] = 1.0
    return mask, pathway_names, subsetted


def sample_negative_binomial(mu: np.ndarray, sigma_per_gene: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """mu: spots x genes. sigma_per_gene: genes, (gamlss NBI dispersion: Var = mu + sigma*mu^2).
    Converts to numpy's negative_binomial(n, p) parameterization: n = 1/sigma (size),
    p = n / (n + mu) (success probability), mean = n(1-p)/p = mu, Var = mu + mu^2/n = mu + sigma*mu^2.
    """
    sigma_per_gene = np.clip(sigma_per_gene, 1e-6, None)
    n = 1.0 / sigma_per_gene  # genes
    n_bcast = np.broadcast_to(n, mu.shape)
    p = n_bcast / (n_bcast + mu)
    p = np.clip(p, 1e-8, 1 - 1e-8)
    return rng.negative_binomial(n_bcast, p)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--backbone-dir", required=True,
                         help="dir with realistic_mu.csv / realistic_sigma.csv from fit_realistic_backbone.R")
    parser.add_argument("--template-dir", required=True,
                         help="dir with the st_coords.csv / st_groups.csv / pathways.gmt the backbone was "
                              "evaluated at (e.g. sim1) -- reused so the injected effect lines up spatially")
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--min-pathway-genes", type=int, default=5,
                         help="drop pathways with fewer surviving member genes than this after "
                              "restricting to genes present in the real lymph node data")
    parser.add_argument("--n-side", type=int, default=30)
    parser.add_argument("--frac-de", type=float, default=0.2)
    parser.add_argument("--de-effect-size", type=float, default=2.0)
    parser.add_argument("--spatial-smoothness", type=float, default=2.0)
    parser.add_argument("--effect-scale", type=float, default=0.3,
                         help="log-mean multiplier for the injected pathway-activity signal; "
                              "kept modest since it stacks additively per gene across all pathways "
                              "the gene belongs to")
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()
    rng = np.random.default_rng(args.seed)

    mu = pd.read_csv(os.path.join(args.backbone_dir, "realistic_mu.csv"), index_col=0)
    sigma = pd.read_csv(os.path.join(args.backbone_dir, "realistic_sigma.csv"), index_col=0)
    coords = pd.read_csv(os.path.join(args.template_dir, "st_coords.csv"), index_col=0)
    groups_df = pd.read_csv(os.path.join(args.template_dir, "st_groups.csv"), index_col=0)

    spot_ids = mu.index.tolist()
    assert spot_ids == coords.index.tolist() == groups_df.index.tolist(), \
        "backbone spots must exactly match the template grid's spot order"
    genes = mu.columns.tolist()
    groups = groups_df["group"].values
    coords_arr = coords[["x", "y"]].values.astype(float)

    mask, pathway_names, subsetted_gmt = subset_pathway_mask(
        os.path.join(args.template_dir, "pathways.gmt"), genes, args.min_pathway_genes
    )
    print(f"Pathway panel after restricting to {len(genes)} real-data genes: "
          f"{len(pathway_names)} / {len(read_gmt(os.path.join(args.template_dir, 'pathways.gmt')))} pathways kept "
          f"(>= {args.min_pathway_genes} surviving member genes)")

    true_activity, de_pathways = simulate_true_activity(
        pathway_names, coords_arr, groups, args.n_side,
        args.frac_de, args.de_effect_size, args.spatial_smoothness, args.seed,
    )
    assert true_activity.index.tolist() == spot_ids

    weights = rng.uniform(0.5, 1.5, size=mask.shape) * mask  # genes x pathways
    log_effect = args.effect_scale * (true_activity.values @ weights.T)  # spots x genes
    mu_final = mu.values * np.exp(log_effect)

    counts = sample_negative_binomial(mu_final, sigma.values.mean(axis=0), rng)  # sigma constant per gene
    counts_df = pd.DataFrame(counts, index=spot_ids, columns=genes)
    expression_df = np.log1p(counts_df)

    os.makedirs(args.out_dir, exist_ok=True)
    expression_df.to_csv(os.path.join(args.out_dir, "st_expression.csv"))
    counts_df.to_csv(os.path.join(args.out_dir, "st_counts_raw.csv"))
    coords.to_csv(os.path.join(args.out_dir, "st_coords.csv"))
    groups_df.to_csv(os.path.join(args.out_dir, "st_groups.csv"))
    true_activity.to_csv(os.path.join(args.out_dir, "true_activity.csv"))
    pd.Series(de_pathways, name="pathway").to_csv(
        os.path.join(args.out_dir, "true_de_pathways.csv"), index=False)
    write_gmt(subsetted_gmt, os.path.join(args.out_dir, "pathways.gmt"))

    print(f"Simulated {len(spot_ids)} spots x {len(genes)} genes, {len(pathway_names)} pathways "
          f"({len(de_pathways)} truly DE between groups) on the real-lymph-node-grounded NB backbone. "
          f"Mean count: {counts.mean():.2f}, frac zero: {(counts == 0).mean():.3f}. "
          f"Written to {args.out_dir}")


if __name__ == "__main__":
    main()
