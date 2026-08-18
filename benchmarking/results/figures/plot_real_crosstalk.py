"""Publication figure for the real-data crosstalk estimate (estimate_real_crosstalk.py):
does the germinal center (GC) micro-environment in real 10x human lymph node Visium data
show more multi-pathway crosstalk than a same-sized random piece of the same tissue?

Not a categorical-identity comparison (this isn't "GRAPHIST vs STAN"), so this figure uses
a standard null-distribution convention instead: a neutral gray histogram (the reference/
background distribution) plus one accent-colored, directly-labeled line for the single
real observation -- satisfies the dataviz skill's "identity never color-alone" principle
via the always-visible annotation rather than a validated categorical pair (the categorical
gate is scoped to identity-bearing series, which this single-observation-vs-null-histogram
pattern isn't).

Outputs: real_crosstalk.pdf (vector) and real_crosstalk.png (300 dpi).
"""
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = os.path.dirname(__file__)
RESULTS_DIR = os.path.join(HERE, "..", "real_crosstalk")

ACCENT = "#2a78d6"  # same blue as GRAPHIST in the dose-response figure -- this analysis
                    # motivates that figure's premise, not a competing data series
NULL_GRAY = "#8a8a86"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 10.5,
    "axes.edgecolor": "#333333",
    "axes.linewidth": 0.8,
    "axes.labelcolor": "#111111",
    "text.color": "#111111",
    "xtick.color": "#333333",
    "ytick.color": "#333333",
    "pdf.fonttype": 42,
})


def clean_axes(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.set_yticks([])
    ax.tick_params(length=3, width=0.8)


def panel(ax, null_vals, real_val, z, pctile, xlabel, title, fmt="{:.3f}"):
    ax.hist(null_vals, bins=18, color=NULL_GRAY, alpha=0.55, edgecolor="white", linewidth=0.6,
             label="random same-sized non-GC subsets (n=200)")
    ax.axvline(real_val, color=ACCENT, linewidth=2.4, zorder=5)
    # Pad the axis so the annotation (always placed left of the line, since real_val sits
    # well right of the null range in this dataset) has room and never clips the figure edge.
    span = real_val - null_vals.min()
    ax.set_xlim(null_vals.min() - span * 0.05, real_val + span * 0.28)
    ymax = ax.get_ylim()[1]
    ax.annotate(f"real GC spots\n{fmt.format(real_val)}\n(z = {z:.1f}, {pctile:.0f}th pct.)",
                 xy=(real_val, ymax * 0.97), xytext=(real_val - span * 0.02, ymax * 0.97),
                 ha="right", va="top", fontsize=8.8, color=ACCENT, fontweight="bold",
                 xycoords="data")
    clean_axes(ax)
    ax.set_xlabel(xlabel)
    ax.set_title(title, loc="left", fontsize=11, fontweight="bold")


def main():
    null_df = pd.read_csv(os.path.join(RESULTS_DIR, "null_distribution.csv"))
    real = pd.read_csv(os.path.join(RESULTS_DIR, "real_gc_stat.csv"), index_col=0).iloc[:, 0]

    real_mean_r = float(real["mean_abs_r"])
    real_n_sig = float(real["n_sig_pairs"])
    z_r = (real_mean_r - null_df["mean_abs_r"].mean()) / null_df["mean_abs_r"].std()
    pct_r = (null_df["mean_abs_r"] < real_mean_r).mean() * 100
    z_n = (real_n_sig - null_df["n_sig_pairs"].mean()) / null_df["n_sig_pairs"].std()
    pct_n = (null_df["n_sig_pairs"] < real_n_sig).mean() * 100

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(9.6, 3.9), constrained_layout=True)

    panel(axA, null_df["mean_abs_r"].values, real_mean_r, z_r, pct_r,
          "Mean |pathway-pathway correlation|\n(disjoint-gene pairs, GSVA scores)",
          "A   Pathway crosstalk strength")
    panel(axB, null_df["n_sig_pairs"].values, real_n_sig, z_n, pct_n,
          "Number of significant pathway pairs\n(BH-corrected, α = 0.05)",
          "B   Pathway crosstalk breadth", fmt="{:.0f}")

    fig.suptitle("Real germinal-center spots show significantly more multi-pathway crosstalk\n"
                  "than a random same-sized region of the same lymph node tissue",
                  fontsize=10.5, y=1.08, color="#111111")

    fig.savefig(os.path.join(HERE, "real_crosstalk.pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(HERE, "real_crosstalk.png"), dpi=300, bbox_inches="tight")
    print("Wrote real_crosstalk.pdf and real_crosstalk.png")


if __name__ == "__main__":
    main()
