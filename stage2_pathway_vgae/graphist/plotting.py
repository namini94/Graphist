"""Plotting helpers for pathway-activity results.

``volcano`` is ported near-verbatim from ``STPA.py`` (already generic/parameterized).
``plot_pathway_spatial`` and ``plot_pathway_boxplot`` generalize the copy-pasted,
hardcoded-pathway-list plotting blocks that appeared at the bottom of each of the
four original scripts into reusable functions driven by a list of pathway names
supplied via dataset config, instead of being hardcoded per file.
"""
from typing import Dict, List, Optional, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from anndata import AnnData
from matplotlib import rcParams


def volcano(
    dfe_res: pd.DataFrame,
    sig_lvl: float = 3.0,
    metric_lvl: float = 3.0,
    annotate_pathways: Optional[Union[str, list]] = None,
    s: int = 10,
    fontsize: int = 10,
    textsize: int = 8,
    figsize=None,
    title: Optional[str] = None,
    save: Optional[str] = None,
):
    """Volcano plot of differential pathway-activity results (Bayes factor vs. mean absolute difference).

    Run :meth:`~graphist.trainer.GraphistTrainer.pathway_differential_activity`
    first to produce ``dfe_res``.

    Parameters
    ----------
    dfe_res
        differential-activity result DataFrame (must have ``bayes_factor`` and
        ``differential_metric`` columns, indexed by pathway name).
    sig_lvl
        absolute Bayes factor cutoff (>= 0).
    metric_lvl
        mean absolute difference cutoff (>= 0).
    annotate_pathways
        pathway name(s) to label. If None, all pathways passing both
        significance thresholds are labeled.
    save
        path to save the figure to (format inferred from extension).
    """
    mad = np.abs(dfe_res["differential_metric"])
    xlim_v = np.abs(dfe_res["bayes_factor"]).max() + 0.5
    ylim_v = mad.max() + 0.5

    idx_sig = np.arange(len(dfe_res["bayes_factor"]))[
        (np.abs(dfe_res["bayes_factor"]) > sig_lvl) & (mad > metric_lvl)
    ]

    fig, ax = plt.subplots(figsize=figsize)
    ax.scatter(dfe_res["bayes_factor"], mad, color="grey", s=s, alpha=0.8, linewidth=0)
    ax.scatter(
        dfe_res["bayes_factor"].values[idx_sig], mad.values[idx_sig], color="red", s=s * 2, linewidth=0
    )
    ax.vlines(x=-sig_lvl, ymin=-0.5, ymax=ylim_v, color="black", linestyles="--", linewidth=1.0, alpha=0.2)
    ax.vlines(x=sig_lvl, ymin=-0.5, ymax=ylim_v, color="black", linestyles="--", linewidth=1.0, alpha=0.2)
    ax.hlines(y=metric_lvl, xmin=-xlim_v, xmax=xlim_v, color="black", linestyles="--", linewidth=1.0, alpha=0.2)

    labels = annotate_pathways if annotate_pathways else dfe_res.index.values[idx_sig]
    if isinstance(labels, str):
        labels = [labels]
    for name in labels:
        i = list(dfe_res.index.values).index(name)
        ax.text(
            x=dfe_res["bayes_factor"].iloc[i],
            y=mad.iloc[i],
            s=name,
            fontdict={"size": textsize},
        )

    ax.set_xlabel(r"$\log_e$(Bayes factor)", fontsize=fontsize)
    ax.set_ylabel("|Differential Metric|", fontsize=fontsize)
    ax.set_ylim([0, ylim_v])
    ax.set_xlim([-xlim_v, xlim_v])
    if title:
        ax.set_title(title, fontsize=fontsize)
    ax.tick_params(axis="x", labelsize=fontsize)
    ax.tick_params(axis="y", labelsize=fontsize)
    plt.grid(False)
    if save:
        plt.savefig(save, format=save.split(".")[-1], dpi=rcParams["savefig.dpi"], bbox_inches="tight")
    plt.show()


def plot_pathway_spatial(
    adata: AnnData,
    pathway_encoded_df: pd.DataFrame,
    pathway: str,
    basis: str = "spatial",
    cmap: str = "coolwarm",
    **kwargs,
):
    """Plot one pathway's activity spatially over the tissue.

    ``basis="spatial"`` uses ``sc.pl.spatial`` (Visium datasets with an H&E image);
    ``basis="coord"`` uses ``sc.pl.embedding`` (datasets without an image, e.g. PDAC).
    """
    adata = adata.copy()
    adata.obs[pathway] = pathway_encoded_df[pathway].values
    if basis == "spatial":
        sc.pl.spatial(adata, color=pathway, cmap=cmap, **kwargs)
    else:
        sc.pl.embedding(adata, basis=basis, color=pathway, cmap=cmap, **kwargs)


def plot_pathway_boxplot(
    adata: AnnData,
    pathway_encoded_df: pd.DataFrame,
    pathways: List[str],
    groupby: str,
    palette: Optional[Dict[str, str]] = None,
    title: str = "Pathway Activity Across Annotations",
    figsize=(10, 6),
):
    """Boxplot of pathway activity scores, grouped by an annotation column (e.g. tissue region)."""
    df = pd.DataFrame(index=adata.obs_names)
    for p in pathways:
        df[p] = pathway_encoded_df[p].values
    df[groupby] = adata.obs[groupby].values
    df_melted = df.melt(id_vars=groupby, var_name="Pathway", value_name="Activity")

    plt.figure(figsize=figsize)
    sns.boxplot(
        data=df_melted,
        x="Activity",
        y="Pathway",
        hue=groupby,
        palette=palette,
        dodge=True,
        fliersize=0.5,
    )
    plt.title(title)
    plt.tight_layout()
    plt.show()
