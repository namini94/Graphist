"""STAN baseline for Task B: spatially-weighted ridge regression for pathway activity.

STAN was designed for TF activity (D = gene-TF prior matrix); we substitute our
gene-pathway mask as D, which is exactly the substitution STAN's own code supports
(`adata.varm['gene_pw']` is its documented fallback name for a non-TF gene set matrix).
Sourced directly from the cloned repo's `stan` package (pure numpy/pandas/scipy, no
extra dependencies beyond what's already installed) rather than requiring the full
package's H&E-image feature-extraction path, which doesn't apply to synthetic data with
no histology image -- the spatial kernel is built from coordinates alone (STAN supports
this natively; the image-feature kernel is documented as an *optional* enhancement).
"""
import argparse
import importlib.util
import os
import sys

import numpy as np
import pandas as pd
from anndata import AnnData

STAN_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "data", "repos", "STAN")


def _load_module(name: str, path: str):
    """Load a single stan/*.py file directly, bypassing stan/__init__.py (which imports
    bigan.py for H&E image feature extraction -- needs torchvision, not installed, and
    not relevant to synthetic data with no histology image anyway).
    """
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


make_kernel_mod = _load_module("stan_make_kernel", os.path.join(STAN_DIR, "stan", "make_kernel.py"))
make_model_mod = _load_module("stan_make_model", os.path.join(STAN_DIR, "stan", "make_model.py"))
make_kernel = make_kernel_mod.make_kernel
Stan = make_model_mod.Stan
assign_folds = make_model_mod.assign_folds


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", required=True)
    parser.add_argument("--out-csv", required=True)
    parser.add_argument("--lam1", type=float, default=5.0)
    parser.add_argument("--lam2", type=float, default=1.0)
    parser.add_argument("--n-folds", type=int, default=2)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()

    expr = pd.read_csv(os.path.join(args.data_dir, "st_expression.csv"), index_col=0)
    coords = pd.read_csv(os.path.join(args.data_dir, "st_coords.csv"), index_col=0).loc[expr.index]

    import sys as _sys
    _sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "..", "stage2_pathway_vgae"))
    from graphist.data.pathway_mask import read_gmt  # noqa: E402

    pathway_dict = read_gmt(os.path.join(args.data_dir, "pathways.gmt"))
    genes = expr.columns.tolist()
    gene_pw = pd.DataFrame(0.0, index=genes, columns=list(pathway_dict.keys()))
    for pathway, members in pathway_dict.items():
        present = [g for g in members if g in gene_pw.index]
        gene_pw.loc[present, pathway] = 1.0

    adata = AnnData(expr.values.astype(np.float32))
    adata.obs_names = expr.index
    adata.var_names = expr.columns
    adata.obsm["spatial"] = coords.values
    adata.varm["gene_pw"] = gene_pw
    adata.layers["dca"] = expr.values.astype(np.float32)

    n_components = min(100, adata.n_obs - 1)
    make_kernel(adata, X=coords.values, n=n_components, kernel_name="kernel")
    assign_folds(adata, n_folds=args.n_folds, random_seed=args.seed)

    model = Stan(adata, kernel_name="kernel", layer="dca")
    model.fit(fixed_params=dict(lam1=args.lam1, lam2=args.lam2))

    # W_concat: pathways x spots -> spots x pathways, matching evaluate.py's convention
    pathway_activity = pd.DataFrame(model.W_concat.T, index=adata.obs_names, columns=gene_pw.columns)
    os.makedirs(os.path.dirname(args.out_csv), exist_ok=True)
    pathway_activity.to_csv(args.out_csv)
    print(f"Wrote STAN pathway activity predictions {pathway_activity.shape} to {args.out_csv}")


if __name__ == "__main__":
    main()
