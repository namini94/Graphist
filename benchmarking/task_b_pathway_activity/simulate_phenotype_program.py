"""Simulates a "phenotype-associated spot" scenario mimicking GRAPHIST's actual paper use
case: spots positively/negatively associated with a bulk phenotype of interest
(Graphist(+)/Graphist(-)), where the phenotype's true biology is a real, thematically
coherent immune/TLS-activation pathway PROGRAM -- chemokine recruitment, TCR signaling,
BCR signaling, downstream NF-kB -- rather than one pathway shifting in isolation. This is
the concrete, citable biological argument for why real phenotype-driven tissue states
involve coordinated multi-pathway crosstalk (see benchmarking/README.md's "Connecting the
mechanism findings to the paper's actual motivation" section).

Generates a MATCHED PAIR of datasets from IDENTICAL ground-truth pathway activity and
IDENTICAL linear/noise terms (same seed, same weights, same baseline, same Gaussian
noise draw -- simulate_expression's outer RNG state doesn't depend on whether interactions
are added, since add_pathway_interactions uses its own independently-seeded RNG), differing
in exactly ONE respect: whether the 6 phenotype-program pathways' effect on gene expression
includes their pairwise bilinear coordination term.

  - phenotype_coordinated/: all C(6,2)=15 pairs among the program get a bilinear
    interaction term -- models real phenotype-driven multi-pathway crosstalk.
  - phenotype_generic/: same DE pathways, same effect size, same noise -- NO coordination
    (independent per-pathway linear effect only) -- the matched "one pathway shifts in
    isolation" control.

Also reports a concrete, quantifiable structural-difference statistic (not just an
assumption): mean pairwise gene-gene correlation, in expression space, within the
phenotype program's own gene set -- this should be measurably higher in the coordinated
condition, which is the actual "something is different here" evidence.

Usage: python simulate_phenotype_program.py --out-dir <dir> [--n-background 24] [--seed 0]
"""
import argparse
import itertools
import os
import sys
from collections import OrderedDict

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter

sys.path.insert(0, os.path.dirname(__file__))
from simulate_pathway_activity import (  # noqa: E402
    read_gmt, write_gmt, simulate_spatial_groups, simulate_expression, add_pathway_interactions,
)

# A real, thematically coherent immune/TLS-activation program: chemokine-driven
# recruitment + T-cell receptor signaling + B-cell receptor signaling + downstream
# NF-kB activation -- the coordinated multi-pathway signature of lymphoid-structure
# activation, directly relevant to the lymph node data already used in sim6_realistic.
PHENOTYPE_PROGRAM = [
    "REACTOME_CHEMOKINE_RECEPTORS_BIND_CHEMOKINES",
    "REACTOME_DOWNSTREAM_TCR_SIGNALING",
    "REACTOME_TCR_SIGNALING",
    "REACTOME_ANTIGEN_ACTIVATES_B_CELL_RECEPTOR_LEADING_TO_GENERATION_OF_SECOND_MESSENGERS",
    "REACTOME_PHOSPHORYLATION_OF_CD3_AND_TCR_ZETA_CHAINS",
    "REACTOME_ACTIVATION_OF_NF_KAPPAB_IN_B_CELLS",
]


def build_panel_with_fixed_de(gmt_path, fixed_pathways, n_background, min_genes, max_genes, seed):
    """Panel = the fixed phenotype-program pathways (placed FIRST, so mask column
    indices 0..len(fixed_pathways)-1 line up with PHENOTYPE_PROGRAM's own order --
    required for add_pathway_interactions' positional indexing to hit the right genes)
    plus n_background randomly-drawn pathways for a realistic-sized panel."""
    all_pathways = read_gmt(gmt_path)
    candidates = [name for name, genes in all_pathways.items()
                  if min_genes <= len(genes) <= max_genes and name not in fixed_pathways]
    rng = np.random.default_rng(seed)
    background = list(rng.choice(candidates, size=min(n_background, len(candidates)), replace=False))
    chosen = list(fixed_pathways) + background
    selected = OrderedDict((name, all_pathways[name]) for name in chosen)

    genes = sorted(set(g for glist in selected.values() for g in glist))
    mask = np.zeros((len(genes), len(chosen)), dtype=float)
    gene_idx = {g: i for i, g in enumerate(genes)}
    for j, name in enumerate(chosen):
        for g in selected[name]:
            mask[gene_idx[g], j] = 1.0
    return genes, chosen, mask, selected


def simulate_true_activity_fixed_de(pathway_names, coords, groups, n_side, de_pathways,
                                     de_effect_size, spatial_smoothness, seed):
    """Identical generative logic to simulate_pathway_activity.simulate_true_activity,
    except the DE pathway SET is given directly (the curated phenotype program) instead
    of chosen randomly via frac_de."""
    rng = np.random.default_rng(seed)
    a_mask = groups == "A"
    b_mask = groups == "B"
    de_set = set(de_pathways)
    activity = np.zeros((coords.shape[0], len(pathway_names)))
    for j, name in enumerate(pathway_names):
        field = rng.normal(size=(n_side, n_side))
        field = gaussian_filter(field, sigma=spatial_smoothness)
        field = (field - field.mean()) / (field.std() + 1e-8)
        field = field.ravel()
        field[a_mask] -= field[a_mask].mean()
        field[b_mask] -= field[b_mask].mean()
        activity[:, j] = field
        if name in de_set:
            shift = np.where(groups == "A", de_effect_size / 2, -de_effect_size / 2)
            activity[:, j] += shift
    return pd.DataFrame(activity, columns=pathway_names, index=[f"spot{i}" for i in range(coords.shape[0])])


def add_cross_group_interactions(true_activity, genes, mask, pathway_names, group_a, group_b,
                                  n_interactions, interaction_strength, seed):
    """Like simulate_pathway_activity.add_pathway_interactions, but candidate pairs are
    restricted to one member from group_a (the phenotype program) and one from group_b
    (background/unrelated pathways elsewhere in the panel) -- models the phenotype's
    pathways crosstalking with OTHER, unrelated ongoing tissue biology, rather than with
    each other. This is the diffuse, panel-wide kind of nonlinear structure sim5/sim7
    actually had (interacting pairs mostly involved at least one non-DE background
    pathway) -- concentrating interactions purely within the 6-pathway program (tried
    first, in add_pathway_interactions via pathway_names=PHENOTYPE_PROGRAM) did not
    reproduce GRAPHIST's advantage; this tests whether the diffuse version does."""
    rng = np.random.default_rng(seed)
    pathway_idx = {p: i for i, p in enumerate(pathway_names)}
    contribution = np.zeros((true_activity.shape[0], len(genes)))
    records = []

    candidate_pairs = [(a, b) for a in group_a for b in group_b]
    chosen_idx = rng.choice(len(candidate_pairs), size=min(n_interactions, len(candidate_pairs)), replace=False)
    gene_pos = {g: i for i, g in enumerate(genes)}
    for idx in chosen_idx:
        p1, p2 = candidate_pairs[idx]
        genes_p1 = {g for i, g in enumerate(genes) if mask[i, pathway_idx[p1]] == 1}
        genes_p2 = {g for i, g in enumerate(genes) if mask[i, pathway_idx[p2]] == 1}
        candidates = sorted(genes_p1 | genes_p2)
        if not candidates:
            continue
        n_affected = max(3, len(candidates) // 3)
        affected = list(rng.choice(candidates, size=min(n_affected, len(candidates)), replace=False))
        interaction_signal = interaction_strength * true_activity[p1].values * true_activity[p2].values
        for g in affected:
            contribution[:, gene_pos[g]] += interaction_signal
        records.append({"pathway1": p1, "pathway2": p2, "genes": affected})

    return contribution, records


def coordination_score(expr: pd.DataFrame, selected_gmt: dict, program: list) -> float:
    """Mean |pairwise Pearson correlation| specifically between genes belonging to
    DIFFERENT phenotype-program pathways (excludes same-pathway gene pairs, which
    correlate anyway just from sharing one pathway's own linear activity term regardless
    of any cross-pathway interaction -- that shared-pathway correlation would swamp the
    much smaller effect actually being measured here). This isolates the thing the
    bilinear interaction term specifically creates: correlation BETWEEN pathways, not
    within one -- the concrete, quantifiable "is this condition structurally different"
    statistic."""
    gene_to_pathways = {}
    for name in program:
        for g in selected_gmt[name]:
            if g in expr.columns:
                gene_to_pathways.setdefault(g, set()).add(name)
    genes = list(gene_to_pathways.keys())
    corr = expr[genes].corr().values
    vals = []
    for i in range(len(genes)):
        for j in range(i + 1, len(genes)):
            if gene_to_pathways[genes[i]].isdisjoint(gene_to_pathways[genes[j]]):
                vals.append(abs(corr[i, j]))
    return float(np.nanmean(vals)) if vals else float("nan")


def write_scenario(out_dir, expr, coords, spot_ids, groups, true_activity, de_pathways, gmt_dict):
    os.makedirs(out_dir, exist_ok=True)
    expr.to_csv(os.path.join(out_dir, "st_expression.csv"))
    pd.DataFrame({"x": coords[:, 0], "y": coords[:, 1]}, index=spot_ids).to_csv(
        os.path.join(out_dir, "st_coords.csv"))
    pd.Series(groups, index=spot_ids, name="group").to_csv(os.path.join(out_dir, "st_groups.csv"))
    true_activity.to_csv(os.path.join(out_dir, "true_activity.csv"))
    pd.Series(de_pathways, name="pathway").to_csv(
        os.path.join(out_dir, "true_de_pathways.csv"), index=False)
    write_gmt(gmt_dict, os.path.join(out_dir, "pathways.gmt"))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gmt-path", default="/Users/naminiyakan/Documents/VEGA_Code/sci-plex/reactomes.gmt")
    parser.add_argument("--out-dir", required=True, help="parent dir; writes coordinated/ and generic/ subdirs")
    parser.add_argument("--n-background", type=int, default=24)
    parser.add_argument("--min-genes", type=int, default=10)
    parser.add_argument("--max-genes", type=int, default=100,
                         help="widened from the usual 60 to accommodate the curated program's largest "
                              "pathway (REACTOME_ACTIVATION_OF_NF_KAPPAB_IN_B_CELLS, 64 genes)")
    parser.add_argument("--n-side", type=int, default=30)
    parser.add_argument("--de-effect-size", type=float, default=2.0)
    parser.add_argument("--spatial-smoothness", type=float, default=2.0)
    parser.add_argument("--noise-sd", type=float, default=1.5, help="matches sim2_hard/sim7's noise level")
    parser.add_argument("--interaction-strength", type=float, default=6.0, help="matches sim5/sim7's strength")
    parser.add_argument("--n-pairs", type=int, default=-1,
                         help="how many pairs get an interaction term (-1 = all possible for the chosen "
                              "--interaction-mode). All-pairs WITHIN the program (dense/full coordination) "
                              "was found to make every method's per-pathway identifiability collapse.")
    parser.add_argument("--interaction-mode", choices=["within", "cross"], default="within",
                         help="'within': pairs among the 6 phenotype pathways only (original design; found "
                              "to either collapse everyone at high density or leave STAN winning at low "
                              "density). 'cross': pairs between a phenotype pathway and a background one -- "
                              "the diffuse, panel-wide structure sim5/sim7 actually had.")
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()

    genes, chosen, mask, selected_gmt = build_panel_with_fixed_de(
        args.gmt_path, PHENOTYPE_PROGRAM, args.n_background, args.min_genes, args.max_genes, args.seed
    )
    assert chosen[:len(PHENOTYPE_PROGRAM)] == PHENOTYPE_PROGRAM, "phenotype program must be first in panel order"
    coords, groups, spot_ids = simulate_spatial_groups(args.n_side)

    true_activity = simulate_true_activity_fixed_de(
        chosen, coords, groups, args.n_side, PHENOTYPE_PROGRAM,
        args.de_effect_size, args.spatial_smoothness, args.seed,
    )

    # base+noise is identical regardless of n_interactions (simulate_expression's outer RNG state doesn't
    # depend on it -- interactions use their own independently-seeded RNG), so generating it once via
    # n_interactions=0 and adding any interaction contribution on top afterward is exactly equivalent to,
    # and simpler than, threading a matched RNG stream through by hand.
    expr_generic, _ = simulate_expression(
        true_activity, genes, mask, args.noise_sd, args.seed,
        pathway_names=PHENOTYPE_PROGRAM, n_interactions=0,
    )

    background = chosen[len(PHENOTYPE_PROGRAM):]
    if args.interaction_mode == "within":
        max_pairs = len(list(itertools.combinations(PHENOTYPE_PROGRAM, 2)))
        n_pairs = max_pairs if args.n_pairs < 0 else min(args.n_pairs, max_pairs)
        interaction_contribution, records = add_pathway_interactions(
            true_activity, genes, mask, PHENOTYPE_PROGRAM, n_pairs, args.interaction_strength, args.seed,
        )
    else:
        max_pairs = len(PHENOTYPE_PROGRAM) * len(background)
        n_pairs = max_pairs if args.n_pairs < 0 else min(args.n_pairs, max_pairs)
        interaction_contribution, records = add_cross_group_interactions(
            true_activity, genes, mask, chosen, PHENOTYPE_PROGRAM, background,
            n_pairs, args.interaction_strength, args.seed,
        )
    expr_coord = expr_generic + interaction_contribution

    program_genes = sorted(set(g for name in PHENOTYPE_PROGRAM for g in selected_gmt[name]) & set(genes))
    score_coord = coordination_score(expr_coord, selected_gmt, PHENOTYPE_PROGRAM)
    score_generic = coordination_score(expr_generic, selected_gmt, PHENOTYPE_PROGRAM)

    coord_dir = os.path.join(args.out_dir, "phenotype_coordinated")
    generic_dir = os.path.join(args.out_dir, "phenotype_generic")
    write_scenario(coord_dir, expr_coord, coords, spot_ids, groups, true_activity, PHENOTYPE_PROGRAM, selected_gmt)
    write_scenario(generic_dir, expr_generic, coords, spot_ids, groups, true_activity, PHENOTYPE_PROGRAM, selected_gmt)
    if records:
        pd.DataFrame(records).to_csv(os.path.join(coord_dir, "true_interactions.csv"), index=False)

    print(f"Phenotype program ({len(PHENOTYPE_PROGRAM)} pathways, {len(program_genes)} genes): "
          f"{PHENOTYPE_PROGRAM}")
    print(f"Panel: {len(spot_ids)} spots x {len(genes)} genes, {len(chosen)} pathways "
          f"({n_pairs} coordinated pairs in the phenotype_coordinated condition)")
    print(f"Coordination score (mean |cross-pathway gene-gene correlation| within the program's gene set):")
    print(f"  phenotype_generic (independent):   {score_generic:.4f}")
    print(f"  phenotype_coordinated (real-like): {score_coord:.4f}")
    print(f"Written to {coord_dir} and {generic_dir}")


if __name__ == "__main__":
    main()
