#!/usr/bin/env python
"""Compare a refactored graphist pipeline run against a saved baseline.

Not a pytest test (needs real data + real compute time) -- a manual
regression check. Usage:

    python validation/reproduce.py --config configs/pdac.yaml \\
        --baseline-dir /path/to/baseline/latent \\
        --output-dir /path/to/scratch/output
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from graphist.clustering import adjusted_rand_index, load_mclust_labels  # noqa: E402
from graphist.config import DatasetConfig  # noqa: E402
from graphist.pipeline import run_pipeline  # noqa: E402


def compare_pathway_matrices(refactored_df: pd.DataFrame, baseline_df: pd.DataFrame) -> None:
    print(f"  refactored shape: {refactored_df.shape}, baseline shape: {baseline_df.shape}")
    if refactored_df.shape != baseline_df.shape:
        print("  SHAPE MISMATCH -- skipping per-column correlation")
        return

    # Baseline CSVs have an unnamed index column (spot names); align on that.
    common_cols = [c for c in baseline_df.columns if c in refactored_df.columns]
    print(f"  {len(common_cols)}/{baseline_df.shape[1]} pathway columns match by name")

    corrs = []
    for col in common_cols:
        a = refactored_df[col].values.astype(float)
        b = baseline_df[col].values.astype(float)
        if a.std() < 1e-8 or b.std() < 1e-8:
            continue
        corrs.append(np.corrcoef(a, b)[0, 1])
    corrs = np.array(corrs)
    print(f"  per-pathway Pearson correlation (refactored vs. baseline): "
          f"mean={corrs.mean():.4f}, median={np.median(corrs):.4f}, "
          f"min={corrs.min():.4f}, frac>0.9={np.mean(corrs > 0.9):.2%}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--baseline-dir", required=True, help="Directory with the baseline latent/*.csv files")
    parser.add_argument("--output-dir", required=True, help="Scratch directory for this run's outputs")
    parser.add_argument("--epochs", type=int, default=None)
    args = parser.parse_args()

    config = DatasetConfig.from_yaml(args.config)
    config.analysis.output_dir = args.output_dir
    if args.epochs is not None:
        config.train.epochs = args.epochs

    print(f"=== Running refactored pipeline for {config.name} ({config.train.epochs} epochs) ===")
    result = run_pipeline(config, save_outputs=True)
    adata = result["adata"]
    print(f"Trained on {adata.n_obs} spots x {adata.n_vars} genes.")

    baseline_pathway_df = pd.read_csv(
        os.path.join(args.baseline_dir, "pathway_encoded_df.csv"), index_col=0
    )
    print("\n--- pathway_encoded_df comparison ---")
    compare_pathway_matrices(result["pathway_encoded_df"], baseline_pathway_df)

    mclust_path = os.path.join(args.baseline_dir, "Mclust_labels.csv")
    if os.path.exists(mclust_path) and config.data.annotation_column in adata.obs.columns:
        labels = load_mclust_labels(mclust_path)
        adata.obs["labels"] = pd.Categorical(labels)
        ari = adjusted_rand_index(adata.obs[config.data.annotation_column], adata.obs["labels"])
        print(f"\nARI (Mclust vs. {config.data.annotation_column}), refactored pipeline: {ari:.4f}")


if __name__ == "__main__":
    main()
