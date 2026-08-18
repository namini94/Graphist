"""Pseudo-bulk phenotype simulation on osmFISH / STARmap, following SpaPheno's recipe.

Real spatial data with known cell-type/layer labels is used two ways:
1. A subsample of cells is repeatedly resampled at controlled mixing proportions of a
   "positive" vs "negative" group to build synthetic pseudo-bulk RNA-seq profiles with a
   known, continuous phenotype value (the positive-group fraction).
2. The full, un-resampled dataset (every cell, all groups) is treated as the "ST data":
   its per-cell group membership is the ground truth a Stage-1 method (GRAPHIST, Scissor,
   SpaPheno, ...) is scored against after being given only the pseudo-bulk cohort + phenotype.

Outputs per scenario (consumed by the R baseline/GRAPHIST runners):
  <out_dir>/bulk_expression.csv   -- genes x pseudo-bulk-samples
  <out_dir>/bulk_phenotype.csv    -- one row per pseudo-bulk sample: phenotype value
  <out_dir>/st_expression.csv     -- genes x real cells (all of them, unresampled)
  <out_dir>/st_coords.csv         -- x,y per real cell
  <out_dir>/st_ground_truth.csv   -- per real cell: -1 negative-group / 0 other / 1 positive-group
  <out_dir>/st_omega.csv          -- cells x cells spatial k-NN adjacency (0/1), built with the
                                      SAME function (graphist.data.graph.generate_adj_mat) used
                                      throughout Stage 2, rather than a second, independent KNN
                                      implementation -- this is what GRAPHIST's Stage 1 R runner
                                      reads as its spatial network regularizer (Omega).
"""
import argparse
import os
import sys
from dataclasses import dataclass

sys.path.insert(
    0, os.path.join(os.path.dirname(__file__), "..", "..", "stage2_pathway_vgae")
)
from graphist.data.graph import generate_adj_mat  # noqa: E402
from anndata import AnnData  # noqa: E402
from typing import List, Optional

import h5py
import numpy as np
import pandas as pd


@dataclass
class LoomDataset:
    genes: List[str]
    expr: np.ndarray  # cells x genes
    labels: np.ndarray  # per-cell group label (string) -- defines the pos/neg phenotype groups
    x: np.ndarray
    y: np.ndarray
    cell_ids: List[str]
    celltype: np.ndarray  # per-cell granular cell-type label -- SpaPheno's composition feature space


def _decode(arr):
    return np.array([v.decode() if isinstance(v, bytes) else v for v in arr])


def load_osmfish(path: str, label_col: str = "Region") -> LoomDataset:
    """``label_col``: "Region" (anatomical layer/region) or "ClusterName" (cell type).

    ``celltype`` is always ClusterName (the finest-grained label available), used as
    SpaPheno's composition feature space regardless of which column defines the
    phenotype groups for a given scenario.
    """
    with h5py.File(path, "r") as f:
        genes = list(_decode(f["row_attrs"]["Gene"][:]))
        expr = f["matrix"][:].T.astype(float)  # -> cells x genes
        labels = _decode(f["col_attrs"][label_col][:])
        celltype = _decode(f["col_attrs"]["ClusterName"][:])
        region = _decode(f["col_attrs"]["Region"][:])
        x = f["col_attrs"]["X"][:]
        y = f["col_attrs"]["Y"][:]
        cell_ids = list(_decode(f["col_attrs"]["CellID"][:]))
    keep = region != "Excluded"  # "Excluded" is always defined on Region, not ClusterName
    return LoomDataset(genes, expr[keep], labels[keep], x[keep], y[keep],
                        [c for c, k in zip(cell_ids, keep) if k], celltype[keep])


def load_starmap(path: str, batch: int = 0) -> LoomDataset:
    """No separate coarse/fine label pair is available for STARmap; ``celltype`` falls
    back to the same "Clusters" label used to define phenotype groups.
    """
    with h5py.File(path, "r") as f:
        genes = list(_decode(f["row_attrs"]["Gene"][:]))
        expr = f["matrix"][:].T.astype(float)
        labels = _decode(f["col_attrs"]["Clusters"][:])
        coords = f["col_attrs"]["Spatial_coordinates"][:]
        batch_id = f["col_attrs"]["BatchID"][:].ravel()
    keep = (labels != "NA") & (batch_id == batch)
    cell_ids = [f"cell{i}" for i in range(expr.shape[0])]
    return LoomDataset(genes, expr[keep], labels[keep], coords[keep, 0], coords[keep, 1],
                        [c for c, k in zip(cell_ids, keep) if k], labels[keep])


def build_pseudobulk_cohort(
    ds: LoomDataset,
    pos_values: List[str],
    neg_values: List[str],
    n_bulk_samples: int = 100,
    cells_per_sample: int = 200,
    seed: int = 0,
) -> tuple:
    """Returns (bulk_expr [samples x genes], bulk_composition [samples x celltypes],
    phenotype [samples], ground_truth [cells] in {-1,0,1}).

    ``bulk_composition`` is the cell-type proportion of the exact same resampled cells
    used to build ``bulk_expr`` for that pseudo-bulk sample (SpaPheno's native feature
    space) -- built from the same draw, not a separate simulation, so the gene-expression
    cohort (GRAPHIST/Scissor) and composition cohort (SpaPheno) are perfectly consistent.
    """
    rng = np.random.default_rng(seed)
    pos_idx = np.where(np.isin(ds.labels, pos_values))[0]
    neg_idx = np.where(np.isin(ds.labels, neg_values))[0]
    assert len(pos_idx) > 0 and len(neg_idx) > 0, "pos/neg groups must be non-empty"

    all_celltypes = sorted(set(ds.celltype))

    bulk_rows = []
    composition_rows = []
    phenotypes = []
    for _ in range(n_bulk_samples):
        frac_pos = rng.uniform(0.0, 1.0)
        n_pos = int(round(cells_per_sample * frac_pos))
        n_neg = cells_per_sample - n_pos
        chosen_pos = rng.choice(pos_idx, size=n_pos, replace=True) if n_pos > 0 else np.array([], dtype=int)
        chosen_neg = rng.choice(neg_idx, size=n_neg, replace=True) if n_neg > 0 else np.array([], dtype=int)
        chosen = np.concatenate([chosen_pos, chosen_neg]).astype(int)
        bulk_rows.append(ds.expr[chosen].mean(axis=0))
        phenotypes.append(frac_pos)

        chosen_types = ds.celltype[chosen]
        counts = pd.Series(chosen_types).value_counts()
        composition_rows.append([counts.get(ct, 0) / len(chosen) for ct in all_celltypes])

    bulk_expr = pd.DataFrame(np.array(bulk_rows), columns=ds.genes,
                              index=[f"pseudobulk{i}" for i in range(n_bulk_samples)])
    bulk_composition = pd.DataFrame(np.array(composition_rows), columns=all_celltypes, index=bulk_expr.index)
    phenotype = pd.Series(phenotypes, index=bulk_expr.index, name="positive_fraction")

    ground_truth = np.zeros(len(ds.labels), dtype=int)
    ground_truth[np.isin(ds.labels, pos_values)] = 1
    ground_truth[np.isin(ds.labels, neg_values)] = -1
    ground_truth = pd.Series(ground_truth, index=ds.cell_ids, name="ground_truth")

    return bulk_expr, bulk_composition, phenotype, ground_truth


# Scenarios, ordered easy -> hard, chosen to mirror SpaPheno's own protocol (including their
# hardest case: subtly-different sub-regions of the same cortical layer).
SCENARIOS = {
    "osmfish_easy": dict(loader="osmfish", label_col="Region", pos=["Layer 4"], neg=["Layer 6"]),
    "osmfish_medium": dict(loader="osmfish", label_col="ClusterName",
                            pos=["Pyramidal L2-3", "Pyramidal L2-3 L5"], neg=["Inhibitory Vip"]),
    "osmfish_hard": dict(loader="osmfish", label_col="Region",
                          pos=["Layer 2-3 lateral"], neg=["Layer 2-3 medial"]),
    "starmap_easy": dict(loader="starmap", pos=["eL2/3"], neg=["eL6-1", "eL6-2"]),
    "starmap_hard": dict(loader="starmap", pos=["eL6-1"], neg=["eL6-2"]),
}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--scenario", required=True, choices=list(SCENARIOS.keys()))
    parser.add_argument("--data-dir", default=os.path.join(os.path.dirname(__file__), "..", "data", "task_a"))
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--n-bulk-samples", type=int, default=100)
    parser.add_argument("--cells-per-sample", type=int, default=200)
    parser.add_argument("--seed", type=int, default=0)
    # k=12 matches the k used throughout Stage 2 (graph_construction(adata, 12) in every
    # dataset config) -- we don't have the literal original neighbors.csv generation
    # parameters for BRCA-PACSI's Stage 1 to confirm this was the same k used there.
    parser.add_argument("--k-neighbors", type=int, default=12)
    args = parser.parse_args()

    spec = SCENARIOS[args.scenario]
    if spec["loader"] == "osmfish":
        ds = load_osmfish(os.path.join(args.data_dir, "osmFISH_SScortex.loom"), label_col=spec["label_col"])
    else:
        ds = load_starmap(os.path.join(args.data_dir, "mpfc_starmap.loom"))

    bulk_expr, bulk_composition, phenotype, ground_truth = build_pseudobulk_cohort(
        ds, spec["pos"], spec["neg"],
        n_bulk_samples=args.n_bulk_samples, cells_per_sample=args.cells_per_sample, seed=args.seed,
    )

    os.makedirs(args.out_dir, exist_ok=True)
    bulk_expr.T.to_csv(os.path.join(args.out_dir, "bulk_expression.csv"))  # genes x samples for R convention
    bulk_composition.to_csv(os.path.join(args.out_dir, "bulk_composition.csv"))  # samples x celltypes
    phenotype.to_csv(os.path.join(args.out_dir, "bulk_phenotype.csv"))
    pd.DataFrame(ds.expr, columns=ds.genes, index=ds.cell_ids).T.to_csv(
        os.path.join(args.out_dir, "st_expression.csv"))  # genes x cells
    pd.DataFrame({"x": ds.x, "y": ds.y}, index=ds.cell_ids).to_csv(os.path.join(args.out_dir, "st_coords.csv"))
    pd.Series(ds.celltype, index=ds.cell_ids, name="celltype").to_csv(os.path.join(args.out_dir, "st_celltype.csv"))
    ground_truth.to_csv(os.path.join(args.out_dir, "st_ground_truth.csv"))

    spatial_adata = AnnData(np.zeros((len(ds.cell_ids), 1)))
    spatial_adata.obsm["spatial"] = np.stack([ds.x, ds.y], axis=1)
    omega = generate_adj_mat(spatial_adata, include_self=False, n=args.k_neighbors)
    pd.DataFrame(omega, index=ds.cell_ids, columns=ds.cell_ids).to_csv(os.path.join(args.out_dir, "st_omega.csv"))

    n_pos = int((ground_truth == 1).sum())
    n_neg = int((ground_truth == -1).sum())
    n_bg = int((ground_truth == 0).sum())
    print(f"Scenario {args.scenario}: {len(ds.genes)} genes, {len(ds.cell_ids)} ST cells "
          f"({n_pos} positive-group, {n_neg} negative-group, {n_bg} other), "
          f"{args.n_bulk_samples} pseudo-bulk samples written to {args.out_dir}")


if __name__ == "__main__":
    main()
