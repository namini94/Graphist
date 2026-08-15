"""High-level orchestration: load -> preprocess -> graph -> train -> save.

This is the config-driven replacement for the procedural script code that
used to sit at the bottom of each of the four original pipeline files. Both
``scripts/run_pipeline.py`` (CLI) and ``validation/reproduce.py`` (regression
checks against the pre-refactor saved outputs) call into this module so the
actual pipeline logic exists in exactly one place.
"""
import os
from typing import Optional

import numpy as np
import pandas as pd
import umap

from .config import DatasetConfig
from .data.graph import graph_construction
from .data.loaders import load_manual_counts_dataset, load_visium_dataset, preprocess
from .data.pathway_mask import create_mask, pathway_names
from .model.graphist_model import GraphistModel
from .trainer import GraphistTrainer
from .utils import seed_everything


def load_dataset(config: DatasetConfig):
    """Load + attach annotations for a dataset, per its DataConfig."""
    dc = config.data
    if dc.loading_mode == "visium":
        adata = load_visium_dataset(dc.data_root)
    elif dc.loading_mode == "manual":
        adata = load_manual_counts_dataset(dc.counts_file, dc.meta_file)
    else:
        raise ValueError(f"Unknown loading_mode: {dc.loading_mode!r}")

    if dc.annotation_file:
        annot = pd.read_csv(dc.annotation_file, sep=dc.annotation_sep, index_col=0)
        if dc.annotation_relabel:
            first_col = annot.columns[0]
            annot[first_col] = annot[first_col].replace(dc.annotation_relabel)
        adata.obs = adata.obs.join(annot)
        if dc.annotation_column not in adata.obs.columns and len(annot.columns) == 1:
            adata.obs[dc.annotation_column] = adata.obs[annot.columns[0]]

    if dc.annotation_column_source and dc.annotation_column_source in adata.obs.columns:
        adata.obs[dc.annotation_column] = adata.obs[dc.annotation_column_source]

    if dc.selection_file:
        selected = pd.read_csv(dc.selection_file, sep=dc.selection_sep, index_col=0, header=0)
        first_col = selected.columns[0]
        selected[first_col] = selected[first_col].replace(dc.selection_relabel)
        adata.obs = adata.obs.join(selected)
        adata.obs[dc.selection_column] = adata.obs[first_col]

    return adata


def run_pipeline(config: DatasetConfig, save_outputs: bool = True) -> dict:
    """Run the full Stage-2 pipeline for one dataset and return its key artifacts.

    Returns a dict with: adata, graph_dict, model, trainer, z, pathway_encoded_df, umap_df.
    """
    seed_everything(config.train.seed)

    adata = load_dataset(config)
    adata = preprocess(
        adata,
        min_cells=config.preprocess.min_cells,
        min_counts=config.preprocess.min_counts,
        target_sum=config.preprocess.target_sum,
        n_top_genes=config.preprocess.n_top_genes,
        use_count_layer=config.preprocess.use_count_layer,
        n_pca_components=config.preprocess.n_pca_components,
        pca_seed=config.preprocess.pca_seed,
    )

    graph_dict = graph_construction(adata, n=config.model.k_neighbors)

    features = adata.var_names.tolist()
    mask = create_mask(
        features,
        config.model.gmt_path,
        add_nodes=config.model.add_nodes,
        min_genes=config.model.min_genes,
        max_genes=config.model.max_genes,
    )
    latent_names = pathway_names(
        config.model.gmt_path,
        add_nodes=config.model.add_nodes,
        min_genes=config.model.min_genes,
        max_genes=config.model.max_genes,
    )
    n_gmvs = len(latent_names)

    model = GraphistModel(
        input_dim=adata.shape[1],
        n_gmvs=n_gmvs,
        mask=mask,
        dropout=config.model.dropout,
        positive_decoder=config.model.positive_decoder,
        gcn_hidden1=config.model.gcn_hidden1,
        p_drop=config.model.p_drop,
    )

    trainer = GraphistTrainer(
        model,
        adata.X,
        graph_dict,
        rec_w=config.train.rec_w,
        gcn_w=config.train.gcn_w,
        gcn_rec_w=config.train.gcn_rec_w,
    )
    trainer.train(epochs=config.train.epochs, lr=config.train.lr, decay=config.train.decay)

    z, de_feat, x_rec = trainer.process()

    reducer = umap.UMAP(random_state=42, min_dist=0.5, n_neighbors=15)
    embedding = reducer.fit_transform(z)
    umap_df = pd.DataFrame(
        {
            "UMAP-1": embedding[:, 0],
            "UMAP-2": embedding[:, 1],
            "annotations": adata.obs[config.data.annotation_column].values
            if config.data.annotation_column in adata.obs.columns
            else None,
        },
        index=adata.obs_names,
    )

    pathway_encoded_df = pd.DataFrame(data=z, index=adata.obs_names, columns=latent_names)

    if save_outputs:
        out_dir = os.path.join(config.analysis.output_dir, "latent")
        os.makedirs(out_dir, exist_ok=True)
        pd.DataFrame(z, index=adata.obs_names).to_csv(os.path.join(out_dir, "z.csv"))
        pathway_encoded_df.to_csv(os.path.join(out_dir, "pathway_encoded_df.csv"))
        umap_df.to_csv(os.path.join(out_dir, "umap.csv"))

    return {
        "adata": adata,
        "graph_dict": graph_dict,
        "model": model,
        "trainer": trainer,
        "z": z,
        "pathway_encoded_df": pathway_encoded_df,
        "umap_df": umap_df,
        "latent_names": latent_names,
    }
