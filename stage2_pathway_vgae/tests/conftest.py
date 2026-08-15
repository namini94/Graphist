"""Shared pytest fixtures: a tiny synthetic spatial dataset + a tiny .gmt file.

Kept deliberately small (40 spots, 60 genes, 5 toy pathways) so the whole
test suite runs in seconds with no GPU and no real data.
"""
import numpy as np
import pandas as pd
import pytest
from anndata import AnnData


N_SPOTS = 40
N_GENES = 60


@pytest.fixture
def rng():
    return np.random.default_rng(0)


@pytest.fixture
def tiny_adata(rng):
    """40 spots on an 8x5 grid, 60 genes, small positive integer counts."""
    xs, ys = np.meshgrid(np.arange(8), np.arange(5))
    coords = np.stack([xs.ravel(), ys.ravel()], axis=1).astype(float)

    counts = rng.poisson(lam=3.0, size=(N_SPOTS, N_GENES)).astype(np.float32)
    var_names = [f"GENE{i}" for i in range(N_GENES)]
    obs_names = [f"spot{i}" for i in range(N_SPOTS)]

    adata = AnnData(
        X=counts,
        obs=pd.DataFrame(index=obs_names),
        var=pd.DataFrame(index=var_names),
    )
    adata.obsm["spatial"] = coords
    adata.obs["group"] = ["A" if i % 2 == 0 else "B" for i in range(N_SPOTS)]
    return adata


@pytest.fixture
def tiny_gmt_path(tmp_path):
    """5 toy pathways covering overlapping subsets of GENE0..GENE59."""
    pathways = {
        "PATHWAY_A": [f"GENE{i}" for i in range(0, 10)],
        "PATHWAY_B": [f"GENE{i}" for i in range(5, 15)],
        "PATHWAY_C": [f"GENE{i}" for i in range(20, 30)],
        "PATHWAY_D": [f"GENE{i}" for i in range(40, 50)],
        "PATHWAY_E": [f"GENE{i}" for i in range(55, 60)],
    }
    path = tmp_path / "toy_pathways.gmt"
    with open(path, "w") as f:
        for name, genes in pathways.items():
            f.write("\t".join([name, "SECOND_COL"] + genes) + "\n")
    return str(path)
