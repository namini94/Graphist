#!/usr/bin/env python
"""CLI entrypoint for the graphist Stage-2 pathway-activity pipeline.

Usage
-----
    python scripts/run_pipeline.py --config configs/pdac.yaml
"""
import argparse
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from graphist.clustering import adjusted_rand_index, load_mclust_labels  # noqa: E402
from graphist.config import DatasetConfig  # noqa: E402
from graphist.pipeline import run_pipeline  # noqa: E402


def main():
    parser = argparse.ArgumentParser(description="Run the graphist Stage-2 pathway-activity pipeline.")
    parser.add_argument("--config", required=True, help="Path to a dataset YAML config (see configs/).")
    parser.add_argument("--epochs", type=int, default=None, help="Override the config's train.epochs.")
    parser.add_argument("--no-save", action="store_true", help="Don't write latent/z.csv etc. to disk.")
    args = parser.parse_args()

    config = DatasetConfig.from_yaml(args.config)
    if args.epochs is not None:
        config.train.epochs = args.epochs

    print(f"Running pipeline for dataset: {config.name}")
    result = run_pipeline(config, save_outputs=not args.no_save)
    adata = result["adata"]
    print(f"Trained on {adata.n_obs} spots x {adata.n_vars} genes -> "
          f"{result['pathway_encoded_df'].shape[1]} pathway-activity latents.")

    if config.analysis.mclust_labels_path:
        labels = load_mclust_labels(config.analysis.mclust_labels_path)
        adata.obs["labels"] = pd.Categorical(labels)
        ari = adjusted_rand_index(adata.obs[config.data.annotation_column], adata.obs["labels"])
        print(f"ARI (Mclust vs. {config.data.annotation_column}): {ari:.4f}")

    if config.analysis.de_enabled:
        gc = config.analysis.de_group_column
        group1_idx = adata.obs_names[adata.obs[gc] == config.analysis.de_group1]
        group2_idx = None
        if config.analysis.de_group2:
            group2_idx = adata.obs_names[adata.obs[gc] == config.analysis.de_group2]

        X = pd.DataFrame(adata.X, columns=adata.var_names, index=adata.obs_names)
        de_df = result["trainer"].pathway_differential_activity(
            X=X,
            group1_idx=group1_idx,
            group2_idx=group2_idx,
            latent_names=result["latent_names"],
            group1_name=str(config.analysis.de_group1),
            group2_name=str(config.analysis.de_group2 or "rest"),
            mode=config.analysis.de_mode,
        )
        out_dir = os.path.join(config.analysis.output_dir, "DEPathways")
        os.makedirs(out_dir, exist_ok=True)
        out_path = os.path.join(out_dir, "differential_activity.csv")
        de_df.to_csv(out_path)
        print(f"Differential pathway-activity results written to {out_path}")
        print(de_df.head(10))


if __name__ == "__main__":
    main()
