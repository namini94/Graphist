"""Publication figure for the Task B dose-response experiment (sim15_dose_*):
sweeps the density of diffuse cross-pathway interaction pairs and shows GRAPHIST's
advantage over STAN emerge, peak, and recede as a function of that single, interpretable
parameter -- the most rigorous (not-a-single-lucky-point) demonstration of the finding in
the whole benchmarking suite. See benchmarking/README.md, "Fifteenth scenario" section.

Palette: dataviz skill's validated default categorical order (fixed slot assignment,
not cycled), first 5 slots -- GRAPHIST/STAN/VEGA/GSVA/ULM in that order:
  blue #2a78d6, orange #eb6834, aqua #1baf7a, yellow #eda100, magenta #e87ba4
Validated via scripts/validate_palette.js: all hard gates pass; the contrast-vs-surface
WARN on 3 slots (aqua/yellow/magenta) is satisfied here via distinct markers/linestyles
and a always-visible legend, per the skill's relief rule (identity is never color-alone).

Outputs (next to this script): dose_response.pdf (vector, for LaTeX/paper submission)
and dose_response.png (300 dpi raster, for slides/preview).
"""
import os

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd

HERE = os.path.dirname(__file__)
DATA_CSV = os.path.join(HERE, "..", "task_b_dose_response.csv")

COLORS = {
    "graphist": "#2a78d6",
    "stan": "#eb6834",
    "vega": "#1baf7a",
    "decoupler_gsva": "#eda100",
    "decoupler_ulm": "#e87ba4",
}
LABELS = {
    "graphist": "GRAPHIST",
    "stan": "STAN",
    "vega": "VEGA (non-spatial ablation)",
    "decoupler_gsva": "GSVA",
    "decoupler_ulm": "ULM",
}
MARKERS = {
    "graphist": "o",
    "stan": "s",
    "vega": "^",
    "decoupler_gsva": "D",
    "decoupler_ulm": "v",
}
# GRAPHIST and STAN are the two methods under direct comparison -- solid lines, full
# weight. The rest are contextual baselines -- lighter weight, dashed, so the eye isn't
# pulled away from the two series the figure is actually about (secondary encoding
# beyond color alone, satisfying the relief rule for the lower-contrast slots).
LINESTYLES = {"graphist": "-", "stan": "-", "vega": "--", "decoupler_gsva": "--", "decoupler_ulm": "--"}
LINEWIDTHS = {"graphist": 2.6, "stan": 2.6, "vega": 1.6, "decoupler_gsva": 1.6, "decoupler_ulm": 1.6}
ZORDER = {"graphist": 5, "stan": 4, "vega": 3, "decoupler_gsva": 2, "decoupler_ulm": 2}
METHOD_ORDER = ["graphist", "stan", "vega", "decoupler_gsva", "decoupler_ulm"]

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
    "svg.fonttype": "none",
    "pdf.fonttype": 42,  # embed as real fonts, not paths -- editable/selectable text in the PDF
})


def clean_axes(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", color="#dddddd", linewidth=0.7, zorder=0)
    ax.set_axisbelow(True)
    ax.tick_params(length=3, width=0.8)


def main():
    df = pd.read_csv(DATA_CSV, index_col=0)
    doses = df.index.values.astype(int)
    # The swept doses (0,5,10,15,25,40,70,100,144) are not evenly spaced, and a linear
    # axis crowds the low end illegibly. Plot at evenly-spaced integer positions and
    # label the ticks with the real dose values instead -- standard fix for an
    # irregularly-sampled sweep; ordering (the only thing that matters here) is preserved.
    x = np.arange(len(doses))

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(9.6, 3.9), constrained_layout=True)

    # ---- Panel A: absolute pathway-recovery accuracy ----
    for m in METHOD_ORDER:
        axA.plot(x, df[m].values, color=COLORS[m], marker=MARKERS[m], markersize=6,
                  linewidth=LINEWIDTHS[m], linestyle=LINESTYLES[m], zorder=ZORDER[m],
                  markeredgecolor="white", markeredgewidth=0.6, label=LABELS[m])
    clean_axes(axA)
    axA.set_xlabel("Diffuse cross-pathway interaction pairs (dose)")
    axA.set_ylabel("Mean Pearson correlation to ground truth")
    axA.set_title("A   Pathway-activity recovery accuracy", loc="left", fontsize=11, fontweight="bold")
    axA.set_ylim(0, 1.0)
    axA.set_xlim(-0.4, len(doses) - 0.6)
    axA.legend(frameon=False, loc="upper right", fontsize=8.8, handlelength=2.2, labelspacing=0.4)

    # ---- Panel B: GRAPHIST's advantage over STAN ----
    gap = df["graphist"] - df["stan"]
    axB.axhline(0, color="#888888", linewidth=0.9, zorder=1)
    axB.fill_between(x, gap.values, 0, where=(gap.values >= 0),
                       color=COLORS["graphist"], alpha=0.18, zorder=1, linewidth=0)
    axB.fill_between(x, gap.values, 0, where=(gap.values < 0),
                       color=COLORS["stan"], alpha=0.18, zorder=1, linewidth=0)
    axB.plot(x, gap.values, color="#333333", linewidth=1.8, zorder=4)
    axB.scatter(x, gap.values, color=np.where(gap.values >= 0, COLORS["graphist"], COLORS["stan"]),
                 s=42, zorder=5, edgecolor="white", linewidth=0.7)

    peak_i = int(np.argmax(gap.values))
    axB.annotate(f"peak advantage\n+{gap.values[peak_i]:.3f}  (n={doses[peak_i]})",
                  xy=(x[peak_i], gap.values[peak_i]), xytext=(x[peak_i] - 1.4, gap.values[peak_i] + 0.035),
                  fontsize=8.6, color="#111111", ha="center",
                  arrowprops=dict(arrowstyle="-", color="#666666", lw=0.8))
    # crossover: first dose where GRAPHIST overtakes STAN
    crossover_i = int(np.argmax(gap.values >= 0))
    axB.annotate(f"crossover\n(n={doses[crossover_i]})", xy=(x[crossover_i], gap.values[crossover_i]),
                  xytext=(x[crossover_i] + 0.3, gap.values[crossover_i] - 0.075),
                  fontsize=8.6, color="#111111", ha="left",
                  arrowprops=dict(arrowstyle="-", color="#666666", lw=0.8))

    axB.text(0.985, 0.94, "favors GRAPHIST", transform=axB.transAxes, ha="right", va="top",
              fontsize=8.6, color=COLORS["graphist"], fontweight="bold")
    axB.text(0.985, 0.07, "favors STAN", transform=axB.transAxes, ha="right", va="bottom",
              fontsize=8.6, color=COLORS["stan"], fontweight="bold")

    clean_axes(axB)
    axB.set_xlabel("Diffuse cross-pathway interaction pairs (dose)")
    axB.set_ylabel("GRAPHIST − STAN  (mean Pearson correlation)")
    axB.set_title("B   GRAPHIST's advantage over STAN", loc="left", fontsize=11, fontweight="bold")
    axB.set_xlim(-0.4, len(doses) - 0.6)
    axB.yaxis.set_major_formatter(mticker.FormatStrFormatter("%+.2f"))

    for ax in (axA, axB):
        ax.set_xticks(x)
        ax.set_xticklabels([str(d) for d in doses])
        ax.tick_params(axis="x", labelrotation=0, labelsize=8.6)

    fig.savefig(os.path.join(HERE, "dose_response.pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(HERE, "dose_response.png"), dpi=300, bbox_inches="tight")
    print("Wrote dose_response.pdf and dose_response.png")


if __name__ == "__main__":
    main()
