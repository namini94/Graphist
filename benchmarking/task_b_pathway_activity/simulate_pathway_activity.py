"""Synthetic ST simulator with known ground-truth per-spot pathway activity.

No real dataset with ground-truth spot-level pathway activity exists in the literature
(see benchmarking/README.md) -- this is the standard workaround: generate expression
FROM known pathway activities using the same gene-pathway structure GRAPHIST's decoder
assumes, so "can a method recover the truth" has an actual answer to check against.

Pipeline:
1. Pick a random subset of real Reactome pathways (from the same reactomes.gmt used by
   the real GRAPHIST pipelines) and their gene universe -- a real, biologically
   structured gene-pathway mask, not an arbitrary synthetic one.
2. Lay out spots on a grid, split into two spatially-contiguous groups (A/B).
3. Sample a spatially-smooth "true activity" surface per pathway (common baseline +
   Gaussian-blurred noise), then inject a known group effect into a designated fraction
   of pathways -- these are the ground-truth differentially-active (DE) pathways.
4. Generate gene expression as a noisy linear combination of pathway activities through
   the gene-pathway mask (the same generative assumption GRAPHIST's masked decoder makes)
   -- this is a "friendly" first test (matches the method's own model class), not an
   adversarial one; harder, non-linear generative variants are a natural follow-up.

Outputs:
  <out_dir>/st_expression.csv    -- spots x genes
  <out_dir>/st_coords.csv        -- x,y per spot
  <out_dir>/st_groups.csv        -- per-spot group label (A/B)
  <out_dir>/true_activity.csv    -- spots x pathways, ground truth
  <out_dir>/true_de_pathways.csv -- list of pathways with an injected group effect
  <out_dir>/pathways.gmt         -- the pathway subset used (so every method sees the
                                     same gene-pathway definitions)
"""
import argparse
import os
from collections import OrderedDict
from typing import List, Tuple

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter


def read_gmt(path: str) -> "OrderedDict[str, list]":
    d = OrderedDict()
    with open(path) as f:
        for line in f:
            parts = line.strip().split("\t")
            d[parts[0]] = parts[2:]
    return d


def write_gmt(d: dict, path: str) -> None:
    with open(path, "w") as f:
        for name, genes in d.items():
            f.write("\t".join([name, "SECOND_COL"] + genes) + "\n")


def build_pathway_mask(
    gmt_path: str, n_pathways: int, min_genes: int, max_genes: int, seed: int
) -> Tuple[List[str], List[str], np.ndarray, dict]:
    """Returns (genes, pathway_names, mask [n_genes x n_pathways], selected_gmt_dict)."""
    all_pathways = read_gmt(gmt_path)
    candidates = [name for name, genes in all_pathways.items() if min_genes <= len(genes) <= max_genes]
    rng = np.random.default_rng(seed)
    chosen = list(rng.choice(candidates, size=min(n_pathways, len(candidates)), replace=False))
    selected = OrderedDict((name, all_pathways[name]) for name in chosen)

    genes = sorted(set(g for glist in selected.values() for g in glist))
    mask = np.zeros((len(genes), len(chosen)), dtype=float)
    gene_idx = {g: i for i, g in enumerate(genes)}
    for j, name in enumerate(chosen):
        for g in selected[name]:
            mask[gene_idx[g], j] = 1.0

    return genes, chosen, mask, selected


def simulate_spatial_groups(n_side: int) -> Tuple[np.ndarray, np.ndarray, List[str]]:
    """n_side x n_side grid; left half = group A, right half = group B."""
    xs, ys = np.meshgrid(np.arange(n_side), np.arange(n_side))
    coords = np.stack([xs.ravel(), ys.ravel()], axis=1).astype(float)
    groups = np.where(coords[:, 0] < n_side / 2, "A", "B")
    spot_ids = [f"spot{i}" for i in range(coords.shape[0])]
    return coords, groups, spot_ids


def simulate_true_activity(
    pathway_names: List[str],
    coords: np.ndarray,
    groups: np.ndarray,
    n_side: int,
    frac_de: float,
    de_effect_size: float,
    spatial_smoothness: float,
    seed: int,
) -> Tuple[pd.DataFrame, List[str]]:
    rng = np.random.default_rng(seed)
    n_pathways = len(pathway_names)
    n_de = max(1, int(round(n_pathways * frac_de)))
    de_pathways = list(rng.choice(pathway_names, size=n_de, replace=False))

    a_mask = groups == "A"
    b_mask = groups == "B"
    activity = np.zeros((coords.shape[0], n_pathways))
    for j, name in enumerate(pathway_names):
        field = rng.normal(size=(n_side, n_side))
        field = gaussian_filter(field, sigma=spatial_smoothness)
        field = (field - field.mean()) / (field.std() + 1e-8)
        field = field.ravel()
        # A spatially-smooth field split into two contiguous halves generically has SOME
        # group-mean difference by chance -- with hundreds of spots per group that becomes
        # "significant" for every pathway, not just the intended DE ones. De-mean within
        # each group separately so every pathway has an EXACT zero group difference by
        # construction; only the designated DE pathways then get a controlled, known shift
        # added back on top. This is what makes true_de_pathways.csv an actual ground truth
        # rather than a fuzzy tendency.
        field[a_mask] -= field[a_mask].mean()
        field[b_mask] -= field[b_mask].mean()
        activity[:, j] = field
        if name in de_pathways:
            shift = np.where(groups == "A", de_effect_size / 2, -de_effect_size / 2)
            activity[:, j] += shift

    df = pd.DataFrame(activity, columns=pathway_names, index=[f"spot{i}" for i in range(coords.shape[0])])
    return df, de_pathways


def simulate_expression(
    true_activity: pd.DataFrame, genes: List[str], mask: np.ndarray, noise_sd: float, seed: int
) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    n_spots = true_activity.shape[0]
    n_genes = len(genes)
    weights = rng.uniform(0.5, 1.5, size=mask.shape) * mask  # gene x pathway, zero outside mask
    baseline = rng.normal(0, 0.3, size=n_genes)
    expr = true_activity.values @ weights.T + baseline
    expr += rng.normal(0, noise_sd, size=expr.shape)
    return pd.DataFrame(expr, columns=genes, index=true_activity.index)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gmt-path", default="/Users/naminiyakan/Documents/VEGA_Code/sci-plex/reactomes.gmt")
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--n-pathways", type=int, default=30)
    parser.add_argument("--min-genes", type=int, default=10)
    parser.add_argument("--max-genes", type=int, default=60)
    parser.add_argument("--n-side", type=int, default=30, help="grid is n_side x n_side spots")
    parser.add_argument("--frac-de", type=float, default=0.2)
    parser.add_argument("--de-effect-size", type=float, default=2.0)
    parser.add_argument("--spatial-smoothness", type=float, default=2.0)
    parser.add_argument("--noise-sd", type=float, default=0.5)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()

    genes, pathway_names, mask, selected_gmt = build_pathway_mask(
        args.gmt_path, args.n_pathways, args.min_genes, args.max_genes, args.seed
    )
    coords, groups, spot_ids = simulate_spatial_groups(args.n_side)
    true_activity, de_pathways = simulate_true_activity(
        pathway_names, coords, groups, args.n_side,
        args.frac_de, args.de_effect_size, args.spatial_smoothness, args.seed,
    )
    expression = simulate_expression(true_activity, genes, mask, args.noise_sd, args.seed)

    os.makedirs(args.out_dir, exist_ok=True)
    expression.to_csv(os.path.join(args.out_dir, "st_expression.csv"))
    pd.DataFrame({"x": coords[:, 0], "y": coords[:, 1]}, index=spot_ids).to_csv(
        os.path.join(args.out_dir, "st_coords.csv"))
    pd.Series(groups, index=spot_ids, name="group").to_csv(os.path.join(args.out_dir, "st_groups.csv"))
    true_activity.to_csv(os.path.join(args.out_dir, "true_activity.csv"))
    pd.Series(de_pathways, name="pathway").to_csv(
        os.path.join(args.out_dir, "true_de_pathways.csv"), index=False)
    write_gmt(selected_gmt, os.path.join(args.out_dir, "pathways.gmt"))

    print(f"Simulated {len(spot_ids)} spots x {len(genes)} genes, {len(pathway_names)} pathways "
          f"({len(de_pathways)} truly DE between groups). Written to {args.out_dir}")


if __name__ == "__main__":
    main()
