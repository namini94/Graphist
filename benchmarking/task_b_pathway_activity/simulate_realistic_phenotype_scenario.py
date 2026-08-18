"""Combines the two things found separately to matter this session:

1. A real-data-grounded noise backbone (fit_realistic_backbone.R -- scDesign3's NB
   marginal model, spatially-smooth mean, fit to real 10x human lymph node data),
   as in sim6_realistic.
2. The specific regime where GRAPHIST beats STAN: a real, curated phenotype-biomarker
   pathway program (chemokine + TCR + BCR + NF-kB signaling) diffusely crosstalking with
   OTHER, unrelated background tissue biology, as in sim14/sim15.

This is the single most biologically-grounded scenario in the whole suite: real
per-gene/per-spot count noise characteristics AND the specific generative structure
that reflects real phenotype-driven biology (a phenotype's biomarker program coupling
to broader tissue pathway activity, not just itself).

Same matched-pair design as simulate_phenotype_program.py: identical ground-truth
activity, identical NB backbone, identical weights -- differing only in whether the
diffuse cross-group interaction term is added before NB sampling.

Requires a backbone already fit on THIS scenario's own gene panel (build_panel_with_fixed_de
with the same --seed/--n-background), since sim6_realistic's backbone was fit on sim1's
panel, which doesn't cover the phenotype program's genes. See fit_realistic_backbone.R,
run against benchmarking/data/task_b/lymph_node_phenotype/{counts,coords}.csv.

Usage: python simulate_realistic_phenotype_scenario.py \
    --backbone-dir <sim16_realistic_phenotype> --template-dir <sim1> --out-dir <out> \
    --n-pairs 40 --interaction-strength 6.0
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from simulate_phenotype_program import (  # noqa: E402
    PHENOTYPE_PROGRAM, build_panel_with_fixed_de, simulate_true_activity_fixed_de,
    add_cross_group_interactions, write_scenario,
)
from simulate_realistic_scenario import sample_negative_binomial  # noqa: E402


def restrict_panel_to_backbone(chosen, mask, selected_gmt, genes, backbone_genes, min_pathway_genes):
    """Restrict the panel to genes present in the real-data-fit backbone (parallels
    simulate_realistic_scenario.subset_pathway_mask), preserving PHENOTYPE_PROGRAM's
    position/order first (required for add_cross_group_interactions' positional
    indexing). Drops pathways left with too few surviving member genes."""
    available = set(backbone_genes)
    kept_pathways, kept_gmt = [], {}
    for name in chosen:
        kept_genes = [g for g in selected_gmt[name] if g in available]
        if len(kept_genes) >= min_pathway_genes:
            kept_pathways.append(name)
            kept_gmt[name] = kept_genes
    new_genes = sorted(set(g for gl in kept_gmt.values() for g in gl))
    gene_idx = {g: i for i, g in enumerate(new_genes)}
    new_mask = np.zeros((len(new_genes), len(kept_pathways)), dtype=float)
    for j, name in enumerate(kept_pathways):
        for g in kept_gmt[name]:
            new_mask[gene_idx[g], j] = 1.0
    return new_genes, kept_pathways, new_mask, kept_gmt


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gmt-path", default="/Users/naminiyakan/Documents/VEGA_Code/sci-plex/reactomes.gmt")
    parser.add_argument("--backbone-dir", required=True)
    parser.add_argument("--template-dir", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--n-background", type=int, default=24)
    parser.add_argument("--min-genes", type=int, default=10)
    parser.add_argument("--max-genes", type=int, default=100)
    parser.add_argument("--min-pathway-genes", type=int, default=5)
    parser.add_argument("--n-side", type=int, default=30)
    parser.add_argument("--de-effect-size", type=float, default=2.0)
    parser.add_argument("--spatial-smoothness", type=float, default=2.0)
    parser.add_argument("--effect-scale", type=float, default=0.3,
                         help="log-mean multiplier for the injected linear pathway-activity signal")
    parser.add_argument("--interaction-effect-scale", type=float, default=0.3,
                         help="log-mean multiplier for the (bounded) interaction contribution")
    parser.add_argument("--interaction-bound", type=float, default=10.0,
                         help="the raw accumulated interaction contribution is heavy-tailed (a phenotype "
                              "pathway can appear in many of the --n-pairs cross-group pairs at once, so "
                              "some genes accumulate contributions from several pairs and can reach ~100+ "
                              "in this simulator's units) -- unlike the additive-Gaussian model sim14/15 "
                              "used, exponentiating an unbounded value in the NB log-link explodes into "
                              "astronomical counts. Squashing with tanh(contribution / bound) * bound caps "
                              "the worst-case log-mean shift at effect_scale * bound regardless of "
                              "accumulation, while preserving relative ranking for typical (small) values.")
    parser.add_argument("--interaction-strength", type=float, default=6.0)
    parser.add_argument("--n-pairs", type=int, default=40,
                         help="defaults to 40 -- the peak-advantage point found in sim15's dose-response sweep")
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()
    rng = np.random.default_rng(args.seed)

    genes, chosen, mask, selected_gmt = build_panel_with_fixed_de(
        args.gmt_path, PHENOTYPE_PROGRAM, args.n_background, args.min_genes, args.max_genes, args.seed
    )
    assert chosen[:len(PHENOTYPE_PROGRAM)] == PHENOTYPE_PROGRAM

    mu = pd.read_csv(os.path.join(args.backbone_dir, "realistic_mu.csv"), index_col=0)
    sigma = pd.read_csv(os.path.join(args.backbone_dir, "realistic_sigma.csv"), index_col=0)
    coords = pd.read_csv(os.path.join(args.template_dir, "st_coords.csv"), index_col=0)
    groups_df = pd.read_csv(os.path.join(args.template_dir, "st_groups.csv"), index_col=0)
    spot_ids = mu.index.tolist()
    assert spot_ids == coords.index.tolist() == groups_df.index.tolist()

    # A handful of genes (9/634 here) get a NaN dispersion from fit_marginal -- almost certainly
    # near-zero-count genes the NB GAMLSS fit couldn't converge a dispersion estimate for. Drop them
    # rather than silently propagating NaN into NB sampling.
    valid_genes = [g for g in mu.columns if not mu[g].isna().any() and not sigma[g].isna().any()]
    if len(valid_genes) < len(mu.columns):
        print(f"Dropping {len(mu.columns) - len(valid_genes)} genes with NaN mu/sigma from fit_marginal: "
              f"{sorted(set(mu.columns) - set(valid_genes))}")
    mu, sigma = mu[valid_genes], sigma[valid_genes]

    genes, chosen, mask, selected_gmt = restrict_panel_to_backbone(
        chosen, mask, selected_gmt, genes, mu.columns.tolist(), args.min_pathway_genes
    )
    phenotype_present = [p for p in PHENOTYPE_PROGRAM if p in chosen]
    assert chosen[:len(phenotype_present)] == phenotype_present, "phenotype program must stay first after restriction"
    mu = mu[genes]
    sigma = sigma[genes]
    groups = groups_df["group"].values
    coords_arr = coords[["x", "y"]].values.astype(float)

    print(f"Phenotype program after restricting to the real backbone's genes: "
          f"{len(phenotype_present)} / {len(PHENOTYPE_PROGRAM)} pathways kept, {len(genes)} genes total, "
          f"{len(chosen)} / {args.n_background + len(PHENOTYPE_PROGRAM)} pathways kept")

    true_activity = simulate_true_activity_fixed_de(
        chosen, coords_arr, groups, args.n_side, phenotype_present,
        args.de_effect_size, args.spatial_smoothness, args.seed,
    )
    assert true_activity.index.tolist() == spot_ids

    weights = rng.uniform(0.5, 1.5, size=mask.shape) * mask  # genes x pathways
    linear_log_effect = args.effect_scale * (true_activity[chosen].values @ weights.T)  # spots x genes

    background = chosen[len(phenotype_present):]
    max_pairs = len(phenotype_present) * len(background)
    n_pairs = min(args.n_pairs, max_pairs)
    interaction_contribution, records = add_cross_group_interactions(
        true_activity, genes, mask, chosen, phenotype_present, background,
        n_pairs, args.interaction_strength, args.seed,
    )

    bounded_interaction = np.tanh(interaction_contribution / args.interaction_bound) * args.interaction_bound
    mu_generic = mu.values * np.exp(linear_log_effect)
    mu_coord = mu.values * np.exp(linear_log_effect + args.interaction_effect_scale * bounded_interaction)

    sigma_per_gene = sigma.values.mean(axis=0)
    counts_generic = sample_negative_binomial(mu_generic, sigma_per_gene, np.random.default_rng(args.seed))
    counts_coord = sample_negative_binomial(mu_coord, sigma_per_gene, np.random.default_rng(args.seed))

    expr_generic = pd.DataFrame(np.log1p(counts_generic), index=spot_ids, columns=genes)
    expr_coord = pd.DataFrame(np.log1p(counts_coord), index=spot_ids, columns=genes)

    coord_dir = os.path.join(args.out_dir, "phenotype_coordinated")
    generic_dir = os.path.join(args.out_dir, "phenotype_generic")
    write_scenario(coord_dir, expr_coord, coords_arr, spot_ids, groups, true_activity, phenotype_present, selected_gmt)
    write_scenario(generic_dir, expr_generic, coords_arr, spot_ids, groups, true_activity, phenotype_present, selected_gmt)
    if records:
        pd.DataFrame(records).to_csv(os.path.join(coord_dir, "true_interactions.csv"), index=False)

    print(f"Simulated {len(spot_ids)} spots x {len(genes)} genes on the real-lymph-node-grounded NB backbone, "
          f"{n_pairs} diffuse cross-group interaction pairs in the coordinated condition.")
    print(f"Mean count -- generic: {counts_generic.mean():.2f} (frac zero {(counts_generic==0).mean():.3f}), "
          f"coordinated: {counts_coord.mean():.2f} (frac zero {(counts_coord==0).mean():.3f})")
    print(f"Written to {coord_dir} and {generic_dir}")


if __name__ == "__main__":
    main()
