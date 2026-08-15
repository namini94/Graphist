"""Dataset loading + preprocessing.

The four original scripts used exactly two loading patterns (manual raw-count
+ metadata text files for PDAC; ``sc.read_visium`` for the other three), and
one shared preprocessing recipe with a handful of per-dataset knobs
(min_cells, whether HVG selection reads from a raw-counts layer or from
already-normalized X, PCA seed/components). Both are captured here so that
per-dataset differences live in config (see ``graphist.config``), not in
copy-pasted code.
"""
from pathlib import Path
from typing import Optional

import pandas as pd
import scanpy as sc
from anndata import AnnData
from sklearn.decomposition import PCA


def load_visium_dataset(data_root: str, make_var_names_unique: bool = True) -> AnnData:
    """Load a 10x Visium sample (the loading pattern used by Maynard/BRCA-COMMOT/BRCA-PACSI)."""
    adata = sc.read_visium(Path(data_root))
    if make_var_names_unique:
        adata.var_names_make_unique()
    return adata


def load_manual_counts_dataset(counts_file: str, meta_file: str, sep: str = "\t") -> AnnData:
    """Load a dataset from separate raw-counts + metadata text files (the PDAC loading pattern).

    ``counts_file`` is expected genes x spots on disk and gets transposed to spots x genes.
    ``meta_file`` must contain ``coor_x``/``coor_y`` columns, which populate both
    ``adata.obsm['coord']`` and ``adata.obsm['spatial']`` (the latter for consistency
    with the Visium loading path so downstream graph construction is loader-agnostic).
    All other metadata columns are attached to ``adata.obs``.
    """
    counts = pd.read_csv(counts_file, sep=sep, index_col=0).T
    meta_df = pd.read_csv(meta_file, sep=sep, index_col=0)

    adata = AnnData(counts)
    coor_df = meta_df.loc[adata.obs_names, ["coor_x", "coor_y"]]
    adata.obsm["coord"] = coor_df.to_numpy()
    adata.obsm["spatial"] = coor_df.to_numpy()
    adata.obs[meta_df.columns] = meta_df.loc[adata.obs_names, meta_df.columns]
    return adata


def preprocess(
    adata: AnnData,
    min_cells: int = 50,
    min_counts: int = 10,
    target_sum: float = 1e6,
    n_top_genes: int = 5000,
    use_count_layer: bool = True,
    n_pca_components: int = 200,
    pca_seed: int = 42,
) -> AnnData:
    """Standard filter -> normalize -> HVG -> scale -> PCA recipe shared by all four datasets.

    ``use_count_layer`` controls whether HVG selection (Seurat v3 flavor, which
    expects raw counts) reads from a ``'count'`` layer snapshotted before
    normalization, or operates directly on already-normalized ``X`` — this
    matches a genuine difference between PDAC (False) and the other three
    datasets (True) in the original scripts, kept here as an explicit,
    documented choice rather than an accidental copy-paste divergence.
    """
    if use_count_layer:
        adata.layers["count"] = adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X.copy()

    sc.pp.filter_genes(adata, min_cells=min_cells)
    sc.pp.filter_genes(adata, min_counts=min_counts)
    sc.pp.normalize_total(adata, target_sum=target_sum)

    hvg_kwargs = {"flavor": "seurat_v3", "n_top_genes": n_top_genes}
    if use_count_layer:
        hvg_kwargs["layer"] = "count"
    sc.pp.highly_variable_genes(adata, **hvg_kwargs)

    adata = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata)

    # sklearn's PCA is used (rather than scanpy's) because scanpy's PCA was found unstable
    # across environments in the original scripts.
    adata.obsm["X_pca"] = PCA(n_components=n_pca_components, random_state=pca_seed).fit_transform(adata.X)

    return adata
