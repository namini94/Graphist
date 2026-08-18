"""Publication figure for the Background-fraction sweep (sim17/sim18): does the
GRAPHIST-vs-STAN advantage survive across a range of realistic-to-unrealistic
Background/+/- tissue compositions, not just the one 71%-Background setting first
tried?

Sweeps patch radius (2,3,4,6,8), which controls the realized Background fraction
(91% down to 36% of the 900-spot grid), at fixed diffuse-crosstalk dose (n_pairs=40,
sim15's peak-advantage setting).
"""
import os

import matplotlib.pyplot as plt
import pandas as pd

HERE = os.path.dirname(__file__)
DATA_CSV = os.path.join(HERE, "..", "task_b_background_sweep.csv")

GRAPHIST_COLOR = "#2a78d6"
STAN_COLOR = "#eb6834"
VEGA_COLOR = "#1baf7a"

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
    ax.grid(axis="y", color="#dddddd", linewidth=0.7, zorder=0)
    ax.set_axisbelow(True)
    ax.tick_params(length=3, width=0.8)


def main():
    df = pd.read_csv(DATA_CSV).sort_values("pct_bg", ascending=False)
    x = df["pct_bg"].values

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(9.6, 3.9), constrained_layout=True)

    axA.plot(x, df["graphist_corr"], "-o", color=GRAPHIST_COLOR, linewidth=2.4, markersize=7,
              markeredgecolor="white", markeredgewidth=0.6, label="GRAPHIST", zorder=4)
    axA.plot(x, df["stan_corr"], "-s", color=STAN_COLOR, linewidth=2.4, markersize=7,
              markeredgecolor="white", markeredgewidth=0.6, label="STAN", zorder=4)
    axA.plot(x, df["vega_corr"], "--^", color=VEGA_COLOR, linewidth=1.6, markersize=6,
              markeredgecolor="white", markeredgewidth=0.6, label="VEGA (ablation)", zorder=3)
    clean_axes(axA)
    axA.set_xlabel("Background fraction of tissue (%)")
    axA.set_ylabel("Full-panel mean Pearson correlation")
    axA.set_title("A   Absolute recovery accuracy", loc="left", fontsize=11, fontweight="bold")
    axA.invert_xaxis()  # realistic (high Background) on the left, reading toward less realistic
    axA.legend(frameon=False, loc="upper left", fontsize=8.8)

    axB.fill_between(x, df["gap_corr"], 0, color=GRAPHIST_COLOR, alpha=0.15, zorder=1)
    axB.axhline(0, color="#888888", linewidth=0.9, zorder=1)
    axB.plot(x, df["gap_corr"], "-o", color="#333333", linewidth=1.8, markersize=7,
              markeredgecolor="white", markeredgewidth=0.6, zorder=4)
    for xi, gi in zip(x, df["gap_corr"]):
        axB.scatter([xi], [gi], color=GRAPHIST_COLOR, s=50, zorder=5, edgecolor="white", linewidth=0.6)
    axB.axvspan(50, x.max() + 5, color=GRAPHIST_COLOR, alpha=0.05, zorder=0)
    axB.text(x.max(), axB.get_ylim()[1] * 0.05, " realistic\n (Background-majority) ", fontsize=7.8,
              color="#666666", ha="left", va="bottom", style="italic")
    clean_axes(axB)
    axB.set_xlabel("Background fraction of tissue (%)")
    axB.set_ylabel("GRAPHIST − STAN  (full-panel correlation)")
    axB.set_title("B   GRAPHIST's advantage over STAN", loc="left", fontsize=11, fontweight="bold")
    axB.invert_xaxis()
    axB.set_ylim(0, max(df["gap_corr"]) * 1.35)

    fig.suptitle("The crosstalk advantage is stable across realistic Background fractions,\n"
                  "eroding only once phenotype-differentiated tissue becomes the majority",
                 x=0.01, ha="left", fontsize=10.8, fontweight="bold", y=1.12)

    fig.savefig(os.path.join(HERE, "background_sweep.pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(HERE, "background_sweep.png"), dpi=300, bbox_inches="tight")
    print("Wrote background_sweep.pdf and background_sweep.png")


if __name__ == "__main__":
    main()
