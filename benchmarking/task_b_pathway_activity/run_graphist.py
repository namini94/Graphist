"""Run GRAPHIST (or its VEGA ablation) on the Task B synthetic dataset.

--mode spatial: the real GRAPHIST pipeline (spatial k-NN graph, GCN encoder, VGAE
    link-reconstruction loss) -- exactly graphist.pipeline's building blocks, just fed
    already-clean synthetic data instead of raw counts (no scanpy preprocessing needed).
--mode nonspatial: the VEGA ablation. GRAPHIST's Stage 2 is explicitly built as a spatial
    extension of VEGA (same pathway-masked decoder, no spatial graph). Rather than writing
    a second model class, this reuses the exact same GraphistModel/GraphistTrainer with an
    identity graph (self-loops only -> the GCN layers degenerate to plain FC layers, since
    spmm(I, X@W) == X@W) and gcn_w=0 (VEGA has no VGAE link-reconstruction term at all,
    since it never had a graph to reconstruct) -- isolating exactly what the spatial graph
    contributes, with everything else (masked decoder, training procedure, loss weights)
    held identical.
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
import torch

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "stage2_pathway_vgae"))
from graphist.data.graph import graph_construction  # noqa: E402
from graphist.data.pathway_mask import create_mask, pathway_names  # noqa: E402
from graphist.model.graphist_model import GraphistModel  # noqa: E402
from graphist.trainer import GraphistTrainer  # noqa: E402
from graphist.utils import seed_everything  # noqa: E402
from anndata import AnnData  # noqa: E402
from sklearn.preprocessing import StandardScaler  # noqa: E402


def identity_graph_dict(n: int) -> dict:
    idx = torch.arange(n)
    indices = torch.stack([idx, idx])
    values = torch.ones(n)
    eye = torch.sparse_coo_tensor(indices, values, (n, n)).coalesce()
    return {"adj_norm": eye, "adj_label": eye, "norm_value": 1.0}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", required=True)
    parser.add_argument("--out-csv", required=True)
    parser.add_argument("--mode", choices=["spatial", "nonspatial"], default="spatial")
    parser.add_argument("--epochs", type=int, default=200)
    parser.add_argument("--k-neighbors", type=int, default=12)
    parser.add_argument("--seed", type=int, default=2024)
    args = parser.parse_args()

    seed_everything(args.seed)

    expr = pd.read_csv(os.path.join(args.data_dir, "st_expression.csv"), index_col=0)
    coords = pd.read_csv(os.path.join(args.data_dir, "st_coords.csv"), index_col=0)
    gmt_path = os.path.join(args.data_dir, "pathways.gmt")

    adata = AnnData(StandardScaler().fit_transform(expr.values))
    adata.obs_names = expr.index
    adata.var_names = expr.columns
    adata.obsm["spatial"] = coords.loc[expr.index].values

    if args.mode == "spatial":
        graph_dict = graph_construction(adata, n=args.k_neighbors)
        gcn_w = 0.1
    else:
        graph_dict = identity_graph_dict(adata.n_obs)
        gcn_w = 0.0  # VEGA has no link-reconstruction term at all

    features = adata.var_names.tolist()
    mask = create_mask(features, gmt_path, add_nodes=1, min_genes=0, max_genes=1000)
    latent_names = pathway_names(gmt_path, add_nodes=1, min_genes=0, max_genes=1000)

    model = GraphistModel(input_dim=adata.n_vars, n_gmvs=len(latent_names), mask=mask,
                           dropout=0.1, positive_decoder=True, gcn_hidden1=800, p_drop=0.2)
    trainer = GraphistTrainer(model, adata.X, graph_dict, rec_w=1000, gcn_w=gcn_w, gcn_rec_w=10)
    trainer.train(epochs=args.epochs, lr=1e-4, decay=5e-4)

    z, _, _ = trainer.process()
    pathway_encoded_df = pd.DataFrame(z, index=adata.obs_names, columns=latent_names)
    os.makedirs(os.path.dirname(args.out_csv), exist_ok=True)
    pathway_encoded_df.to_csv(args.out_csv)
    print(f"[{args.mode}] Wrote pathway activity predictions ({pathway_encoded_df.shape}) to {args.out_csv}")


if __name__ == "__main__":
    main()
