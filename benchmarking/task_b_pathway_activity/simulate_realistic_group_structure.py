"""Tests whether the sim14/15 GRAPHIST-vs-STAN dose-response finding survives a more
realistic spot-group structure.

Every prior Task B scenario (sim1-sim16) modeled the ENTIRE tissue as a clean two-way
split -- group A = one spatially-contiguous half of the grid, group B = the other half,
50/50, nothing else. That does not match what GRAPHIST's actual Stage 1 produces on real
data: a Scissor-style regression typically yields THREE categories -- Graphist(+),
Graphist(-), and Background (spots not significantly associated with the phenotype at
all) -- where Background is often the *majority*, and +/- are smaller, scattered minority
populations, not each occupying one clean contiguous half of the tissue.

This scenario:
1. Lays out three groups on the same 30x30 grid: Background (the majority, no
   phenotype-driven signal at all), and Graphist(+)/Graphist(-) as SCATTERED minority
   patches (several small Gaussian-blob foci each, not one contiguous block) -- closer to
   how real Stage 1 output looks spatially (e.g. several small tumor nests or immune-
   active niches, not half the slide).
2. Injects the known DE shift ONLY into +/- spots for the phenotype-program pathways;
   Background spots get exactly zero shift (the null-phenotype-association baseline).
3. Injects the SAME diffuse cross-group interaction mechanism that gave GRAPHIST its
   dose-response advantage in sim14/15 (at n_pairs=40, the peak-advantage dose found
   there), so the phenotype's biomarker program still crosstalks with background tissue
   biology -- but now that crosstalk, like the DE shift itself, is only really "live" in
   the +/- minority, diluted by a majority Background population with no such structure.
4. Runs Stage 2 (GRAPHIST/VEGA/STAN/decoupleR) on the WHOLE dataset (Background included
   -- matching real usage: Stage 2 infers on everyone, using spatial context from every
   neighboring spot regardless of phenotype status), then differential-tests only
   Graphist(+) vs. Graphist(-) (Background excluded from the DE test itself, matching how
   the real biomarker-discovery step only ever compares + vs -).

Usage: python simulate_realistic_group_structure.py --out-dir <dir> [--n-pairs 40] [--seed 0]
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter

sys.path.insert(0, os.path.dirname(__file__))
from simulate_phenotype_program import (  # noqa: E402
    PHENOTYPE_PROGRAM, build_panel_with_fixed_de, add_cross_group_interactions, write_scenario,
)
from simulate_pathway_activity import simulate_expression, write_gmt  # noqa: E402


def simulate_scattered_groups(n_side: int, n_plus_patches: int, n_minus_patches: int,
                               patch_radius: float, seed: int):
    """Background (majority) + scattered +/- minority patches, instead of one clean
    contiguous half each."""
    rng = np.random.default_rng(seed)
    xs, ys = np.meshgrid(np.arange(n_side), np.arange(n_side))
    coords = np.stack([xs.ravel(), ys.ravel()], axis=1).astype(float)
    groups = np.full(coords.shape[0], "Background", dtype=object)

    def place_patches(n_patches, label, avoid_mask):
        centers = rng.uniform(0, n_side, size=(n_patches, 2))
        for cx, cy in centers:
            dist = np.sqrt((coords[:, 0] - cx) ** 2 + (coords[:, 1] - cy) ** 2)
            hit = (dist < patch_radius) & ~avoid_mask
            groups[hit] = label

    place_patches(n_plus_patches, "+", np.zeros(coords.shape[0], dtype=bool))
    place_patches(n_minus_patches, "-", groups == "+")

    spot_ids = [f"spot{i}" for i in range(coords.shape[0])]
    return coords, groups, spot_ids


def simulate_true_activity_3way(pathway_names, coords, groups, n_side, de_pathways,
                                 de_effect_size, spatial_smoothness, seed):
    """Same smooth-field generative logic as elsewhere, globally de-meaned (not
    per-group, since Background has no assigned 'side' to de-mean against) so every
    pathway has exactly zero built-in group difference by construction; only the
    designated DE pathways get a deterministic +/- shift added for Graphist(+)/(-)
    spots specifically. Background spots get NO shift at all -- the null-association
    baseline, unlike every prior scenario where every spot belonged to a shifted group."""
    rng = np.random.default_rng(seed)
    plus_mask, minus_mask = groups == "+", groups == "-"
    de_set = set(de_pathways)
    activity = np.zeros((coords.shape[0], len(pathway_names)))
    for j, name in enumerate(pathway_names):
        field = rng.normal(size=(n_side, n_side))
        field = gaussian_filter(field, sigma=spatial_smoothness)
        field = (field - field.mean()) / (field.std() + 1e-8)
        field = field.ravel()
        field -= field.mean()  # global de-mean: zero built-in difference for ANY grouping
        activity[:, j] = field
        if name in de_set:
            activity[plus_mask, j] += de_effect_size / 2
            activity[minus_mask, j] -= de_effect_size / 2
            # Background spots: no shift at all -- exactly the null-association baseline
    return pd.DataFrame(activity, columns=pathway_names, index=[f"spot{i}" for i in range(coords.shape[0])])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gmt-path", default="/Users/naminiyakan/Documents/VEGA_Code/sci-plex/reactomes.gmt")
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--n-background", type=int, default=24)
    parser.add_argument("--min-genes", type=int, default=10)
    parser.add_argument("--max-genes", type=int, default=100)
    parser.add_argument("--n-side", type=int, default=30)
    parser.add_argument("--n-plus-patches", type=int, default=4)
    parser.add_argument("--n-minus-patches", type=int, default=4)
    parser.add_argument("--patch-radius", type=float, default=4.0,
                         help="radius (grid units) of each +/- scattered patch")
    parser.add_argument("--de-effect-size", type=float, default=2.0)
    parser.add_argument("--spatial-smoothness", type=float, default=2.0)
    parser.add_argument("--noise-sd", type=float, default=1.5)
    parser.add_argument("--interaction-strength", type=float, default=6.0)
    parser.add_argument("--n-pairs", type=int, default=40,
                         help="diffuse cross-group interaction pairs -- 40 is sim15's peak-advantage dose")
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()

    genes, chosen, mask, selected_gmt = build_panel_with_fixed_de(
        args.gmt_path, PHENOTYPE_PROGRAM, args.n_background, args.min_genes, args.max_genes, args.seed
    )
    coords, groups, spot_ids = simulate_scattered_groups(
        args.n_side, args.n_plus_patches, args.n_minus_patches, args.patch_radius, args.seed
    )
    n_bg, n_plus, n_minus = (groups == "Background").sum(), (groups == "+").sum(), (groups == "-").sum()
    print(f"Spot groups: Background={n_bg}, Graphist(+)={n_plus}, Graphist(-)={n_minus} "
          f"(out of {len(spot_ids)} total)")

    true_activity = simulate_true_activity_3way(
        chosen, coords, groups, args.n_side, PHENOTYPE_PROGRAM,
        args.de_effect_size, args.spatial_smoothness, args.seed,
    )

    expr_base, _ = simulate_expression(
        true_activity, genes, mask, args.noise_sd, args.seed,
        pathway_names=PHENOTYPE_PROGRAM, n_interactions=0,
    )
    background_pathways = chosen[len(PHENOTYPE_PROGRAM):]
    interaction_contribution, records = add_cross_group_interactions(
        true_activity, genes, mask, chosen, PHENOTYPE_PROGRAM, background_pathways,
        args.n_pairs, args.interaction_strength, args.seed,
    )
    expr = expr_base + interaction_contribution

    out_dir = args.out_dir
    os.makedirs(out_dir, exist_ok=True)
    expr.to_csv(os.path.join(out_dir, "st_expression.csv"))
    pd.DataFrame({"x": coords[:, 0], "y": coords[:, 1]}, index=spot_ids).to_csv(
        os.path.join(out_dir, "st_coords.csv"))
    pd.Series(groups, index=spot_ids, name="group").to_csv(os.path.join(out_dir, "st_groups.csv"))
    true_activity.to_csv(os.path.join(out_dir, "true_activity.csv"))
    pd.Series(PHENOTYPE_PROGRAM, name="pathway").to_csv(
        os.path.join(out_dir, "true_de_pathways.csv"), index=False)
    write_gmt(selected_gmt, os.path.join(out_dir, "pathways.gmt"))
    if records:
        pd.DataFrame(records).to_csv(os.path.join(out_dir, "true_interactions.csv"), index=False)

    print(f"Simulated {len(spot_ids)} spots x {len(genes)} genes, {len(chosen)} pathways, "
          f"{args.n_pairs} diffuse crosstalk pairs. Written to {out_dir}")


if __name__ == "__main__":
    main()
