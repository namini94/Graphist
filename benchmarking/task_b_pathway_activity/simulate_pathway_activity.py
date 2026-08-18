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

Two nonlinearity knobs, in increasing order of how fundamentally they break the linear
generative assumption (see --nonlinearity and --n-interactions):
  - "saturating": a monotonic tanh compression. Found empirically to be an insufficient
    stress test -- Spearman correlation stayed higher than Pearson for every method
    tested, meaning rank order survived well enough that linear methods (which a Pearson
    correlation with a monotonically-transformed linear signal still rewards) barely
    noticed. A real result, not a null one: it shows *why* a harder nonlinearity is
    needed, not just that one wasn't tried.
  - interaction terms (--n-interactions): some genes' expression additionally depends on
    the PRODUCT of two pathways' activities, not a linear combination of each
    independently. No fixed linear gene-pathway design matrix D can represent a bilinear
    term Y ~ activity_p1 * activity_p2 for a per-spot-independent linear solve (which is
    exactly how STAN and any other per-spot ridge/regression method operates) --
    structurally, only a method that learns one shared nonlinear function across all
    spots jointly (e.g. GRAPHIST's GCN encoder, whose weights are fit across the whole
    dataset at once) has any architectural path to compensate for it.
"""
import argparse
import itertools
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


def apply_nonlinearity(activity: np.ndarray, mode: str, scale: float) -> np.ndarray:
    """Transform true activity before it enters the linear expression combination.

    ``mode="none"``: identity -- the original "friendly", exactly-linear generative
    process (matches what a masked-linear decoder / linear ridge regression assumes).

    ``mode="saturating"``: ``scale * tanh(activity / scale)`` -- a biologically
    motivated nonlinearity (transcriptional responses to regulatory/pathway signal
    commonly saturate, e.g. Michaelis-Menten-like kinetics), not an adversarial trick.
    Small ``scale`` = strong saturation = a harder deviation from linearity. Breaks the
    exact linear generative assumption that STAN's closed-form ridge regression (and any
    purely linear method) depends on, while ground truth for evaluation purposes stays
    the PRE-saturation activity -- the quantity we actually want a method to recover.
    """
    if mode == "none":
        return activity
    if mode == "saturating":
        return scale * np.tanh(activity / scale)
    raise ValueError(f"Unknown nonlinearity mode: {mode!r}")


def add_pathway_interactions(
    true_activity: pd.DataFrame,
    genes: List[str],
    mask: np.ndarray,
    pathway_names: List[str],
    n_interactions: int,
    interaction_strength: float,
    seed: int,
) -> Tuple[np.ndarray, List[dict]]:
    """Adds a bilinear (product-of-two-pathways) contribution to a subset of genes.

    Returns (contribution [spots x genes], list of {pathway1, pathway2, genes} records
    describing which pairs/genes were affected, for documentation/debugging).

    Structurally different from ``apply_nonlinearity``: that one is still a function of a
    SINGLE pathway's own activity (however nonlinear), so a linear method can still
    partially track it via correlation with that one pathway. A product of TWO
    pathways' activities cannot be written as `sum_p weight_p * activity_p` for any
    choice of per-pathway weights -- no fixed linear design matrix D represents it,
    which is exactly the assumption every method here except GRAPHIST's jointly-trained
    nonlinear encoder is built on.
    """
    rng = np.random.default_rng(seed)
    pathway_idx = {p: i for i, p in enumerate(pathway_names)}
    contribution = np.zeros((true_activity.shape[0], len(genes)))
    records = []

    candidate_pairs = list(itertools.combinations(pathway_names, 2))
    chosen = rng.choice(len(candidate_pairs), size=min(n_interactions, len(candidate_pairs)), replace=False)
    for idx in chosen:
        p1, p2 = candidate_pairs[idx]
        genes_p1 = {g for i, g in enumerate(genes) if mask[i, pathway_idx[p1]] == 1}
        genes_p2 = {g for i, g in enumerate(genes) if mask[i, pathway_idx[p2]] == 1}
        candidates = sorted(genes_p1 | genes_p2)
        if not candidates:
            continue
        n_affected = max(3, len(candidates) // 3)
        affected = list(rng.choice(candidates, size=min(n_affected, len(candidates)), replace=False))

        interaction_signal = interaction_strength * true_activity[p1].values * true_activity[p2].values
        gene_pos = {g: i for i, g in enumerate(genes)}
        for g in affected:
            contribution[:, gene_pos[g]] += interaction_signal
        records.append({"pathway1": p1, "pathway2": p2, "genes": affected})

    return contribution, records


def simulate_expression(
    true_activity: pd.DataFrame,
    genes: List[str],
    mask: np.ndarray,
    noise_sd: float,
    seed: int,
    nonlinearity: str = "none",
    nonlinearity_scale: float = 1.0,
    pathway_names: List[str] = None,
    n_interactions: int = 0,
    interaction_strength: float = 1.0,
) -> Tuple[pd.DataFrame, List[dict]]:
    rng = np.random.default_rng(seed)
    n_spots = true_activity.shape[0]
    n_genes = len(genes)
    weights = rng.uniform(0.5, 1.5, size=mask.shape) * mask  # gene x pathway, zero outside mask
    baseline = rng.normal(0, 0.3, size=n_genes)
    effective_activity = apply_nonlinearity(true_activity.values, nonlinearity, nonlinearity_scale)
    expr = effective_activity @ weights.T + baseline

    interaction_records = []
    if n_interactions > 0:
        interaction_contribution, interaction_records = add_pathway_interactions(
            true_activity, genes, mask, pathway_names, n_interactions, interaction_strength, seed
        )
        expr += interaction_contribution

    expr += rng.normal(0, noise_sd, size=expr.shape)
    return pd.DataFrame(expr, columns=genes, index=true_activity.index), interaction_records


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
    parser.add_argument("--nonlinearity", choices=["none", "saturating"], default="none")
    parser.add_argument("--nonlinearity-scale", type=float, default=1.0)
    parser.add_argument("--n-interactions", type=int, default=0,
                         help="number of pathway pairs given a bilinear (product) interaction term")
    parser.add_argument("--interaction-strength", type=float, default=1.0)
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
    expression, interaction_records = simulate_expression(
        true_activity, genes, mask, args.noise_sd, args.seed,
        nonlinearity=args.nonlinearity, nonlinearity_scale=args.nonlinearity_scale,
        pathway_names=pathway_names, n_interactions=args.n_interactions,
        interaction_strength=args.interaction_strength,
    )

    os.makedirs(args.out_dir, exist_ok=True)
    expression.to_csv(os.path.join(args.out_dir, "st_expression.csv"))
    pd.DataFrame({"x": coords[:, 0], "y": coords[:, 1]}, index=spot_ids).to_csv(
        os.path.join(args.out_dir, "st_coords.csv"))
    pd.Series(groups, index=spot_ids, name="group").to_csv(os.path.join(args.out_dir, "st_groups.csv"))
    true_activity.to_csv(os.path.join(args.out_dir, "true_activity.csv"))
    pd.Series(de_pathways, name="pathway").to_csv(
        os.path.join(args.out_dir, "true_de_pathways.csv"), index=False)
    write_gmt(selected_gmt, os.path.join(args.out_dir, "pathways.gmt"))
    if interaction_records:
        pd.DataFrame(interaction_records).to_csv(os.path.join(args.out_dir, "true_interactions.csv"), index=False)

    print(f"Simulated {len(spot_ids)} spots x {len(genes)} genes, {len(pathway_names)} pathways "
          f"({len(de_pathways)} truly DE between groups, {len(interaction_records)} interacting pairs). "
          f"Written to {args.out_dir}")


if __name__ == "__main__":
    main()
