"""Publication figure for the lymph-node known-biology check: does GRAPHIST's real
GC-vs-non-GC differential-activity call align better with established GC immunology
than STAN's does, on the SAME real tissue?

Grouped bar chart: |Cohen's d| (real GC vs. non-GC) for each of the 6 known GC-relevant
pathways (chemokine recruitment, TCR/BCR signaling, NF-kB activation), GRAPHIST vs. STAN
side by side, with the DE-calling threshold (|d|=0.3) marked -- bars crossing it are
"called DE"; those that don't (STAN's three TCR-related pathways) are the visible gap.
"""
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = os.path.dirname(__file__)
DATA_DIR = os.path.join(HERE, "..", "..", "task_b_pathway_activity")

GRAPHIST_COLOR = "#2a78d6"
STAN_COLOR = "#eb6834"
THRESHOLD = 0.3

# Short display labels for the 6 known GC-immunology pathways
LABELS = {
    "REACTOME_CHEMOKINE_RECEPTORS_BIND_CHEMOKINES": "Chemokine\nreceptors",
    "REACTOME_ANTIGEN_ACTIVATES_B_CELL_RECEPTOR_LEADING_TO_GENERATION_OF_SECOND_MESSENGERS": "BCR\nsignaling",
    "REACTOME_DOWNSTREAM_TCR_SIGNALING": "Downstream\nTCR signaling",
    "REACTOME_TCR_SIGNALING": "TCR\nsignaling",
    "REACTOME_PHOSPHORYLATION_OF_CD3_AND_TCR_ZETA_CHAINS": "CD3/TCR-zeta\nphosphorylation",
    "REACTOME_ACTIVATION_OF_NF_KAPPAB_IN_B_CELLS": "NF-kB activation\n(B cells)",
}

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


def cohens_d_table():
    from scipy.stats import ttest_ind
    groups = pd.read_csv(os.path.join(DATA_DIR, "..", "data", "task_b", "lymph_node_real_gc_check",
                                       "st_groups.csv"), index_col=0)["group"]
    a_idx, b_idx = groups[groups == "A"].index, groups[groups == "B"].index
    out = {}
    for method in ["graphist", "stan"]:
        pred = pd.read_csv(os.path.join(DATA_DIR, "..", "data", "task_b", "lymph_node_real_gc_check",
                                         f"{method}_predictions.csv"), index_col=0)
        d = {}
        for p in LABELS:
            if p not in pred.columns:
                continue
            a, b = pred.loc[a_idx, p].values, pred.loc[b_idx, p].values
            pooled_sd = np.sqrt((a.var(ddof=1) + b.var(ddof=1)) / 2)
            d[p] = (a.mean() - b.mean()) / (pooled_sd + 1e-12)
        out[method] = d
    return out


def main():
    d = cohens_d_table()
    pathways = list(LABELS.keys())
    graphist_vals = np.array([abs(d["graphist"].get(p, 0)) for p in pathways])
    stan_vals = np.array([abs(d["stan"].get(p, 0)) for p in pathways])

    x = np.arange(len(pathways))
    width = 0.36

    fig, ax = plt.subplots(figsize=(8.4, 4.0), constrained_layout=True)
    b1 = ax.bar(x - width / 2, graphist_vals, width, color=GRAPHIST_COLOR, label="GRAPHIST",
                edgecolor="white", linewidth=0.6, zorder=3)
    b2 = ax.bar(x + width / 2, stan_vals, width, color=STAN_COLOR, label="STAN",
                edgecolor="white", linewidth=0.6, zorder=3)

    ax.axhline(THRESHOLD, color="#888888", linewidth=1.0, linestyle="--", zorder=2)
    ax.text(len(pathways) - 0.55, THRESHOLD + 0.03, "DE-calling threshold (|d| = 0.3)",
            fontsize=8.2, color="#666666", ha="right", va="bottom")

    for bars, color, vals in [(b1, GRAPHIST_COLOR, graphist_vals), (b2, STAN_COLOR, stan_vals)]:
        for bar, v in zip(bars, vals):
            called = v > THRESHOLD
            ax.text(bar.get_x() + bar.get_width() / 2, v + 0.04, "DE" if called else "n.s.",
                    ha="center", va="bottom", fontsize=8, color=color if called else "#999999",
                    fontweight="bold" if called else "normal")

    ax.set_xticks(x)
    ax.set_xticklabels([LABELS[p] for p in pathways], fontsize=8.6)
    ax.set_ylabel("|Cohen's d|  (real GC vs. non-GC spots)")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", color="#dddddd", linewidth=0.7, zorder=0)
    ax.set_axisbelow(True)
    ax.set_ylim(0, max(graphist_vals.max(), stan_vals.max()) + 0.3)
    ax.legend(frameon=False, loc="upper right", fontsize=9)
    ax.set_title("GRAPHIST recovers all 6 known GC-immunology pathways as differentially active;\n"
                  "STAN misses all 3 TCR-related ones", loc="left", fontsize=10.8, fontweight="bold")

    fig.savefig(os.path.join(HERE, "gc_known_biology.pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(HERE, "gc_known_biology.png"), dpi=300, bbox_inches="tight")
    print("Wrote gc_known_biology.pdf and gc_known_biology.png")


if __name__ == "__main__":
    main()
