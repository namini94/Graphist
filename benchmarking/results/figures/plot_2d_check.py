"""Confirms the two Task B "dose" axes are orthogonal, not contradictory: the
crosstalk-density dose-response (sim15) and the Background-coverage sweep (sim18) are
different knobs. This overlays the crosstalk-density curve at two different fixed
Background-coverage levels -- the original 50/50 split (0% Background, every spot
phenotype-associated) and the realistic 71%-Background 3-way structure (sim17/19) --
to show the SHAPE (crossover, peak-at-40, decline) is preserved, only the absolute
scale shifts. If the two axes were actually the same thing in disguise, this overlay
would look nothing alike; instead the peak lands at the same dose in both.
"""
import os

import matplotlib.pyplot as plt
import pandas as pd

HERE = os.path.dirname(__file__)

ORIGINAL_GAP = {0: -0.103, 5: 0.006, 10: 0.036, 15: 0.072, 25: 0.103, 40: 0.128,
                70: 0.071, 100: 0.062, 144: 0.056}

GRAPHIST_COLOR = "#2a78d6"

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


def main():
    df = pd.read_csv(os.path.join(HERE, "..", "task_b_2d_check.csv")).sort_values("n_pairs")
    df["gap"] = df["graphist"] - df["stan"]

    orig_x = sorted(ORIGINAL_GAP.keys())
    orig_y = [ORIGINAL_GAP[x] for x in orig_x]

    fig, ax = plt.subplots(figsize=(6.6, 4.2), constrained_layout=True)
    ax.axhline(0, color="#888888", linewidth=0.9, zorder=1)

    ax.plot(orig_x, orig_y, "-o", color="#999999", linewidth=1.8, markersize=6,
            markeredgecolor="white", markeredgewidth=0.6, label="0% Background (original 50/50 split)", zorder=3)
    ax.plot(df["n_pairs"], df["gap"], "-o", color=GRAPHIST_COLOR, linewidth=2.4, markersize=7,
            markeredgecolor="white", markeredgewidth=0.6, label="71% Background (realistic, sim17/19)", zorder=4)

    for x, y in [(40, 0.128), (40, 0.097)]:
        ax.axvline(40, color="#cccccc", linewidth=1.0, linestyle=":", zorder=0)

    ax.set_xlabel("Diffuse cross-pathway interaction pairs (dose)")
    ax.set_ylabel("GRAPHIST − STAN  (mean Pearson correlation)")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", color="#dddddd", linewidth=0.7, zorder=0)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, loc="lower right", fontsize=8.6)
    ax.set_title("Same dose-response shape at two different Background levels:\n"
                  "the two axes are orthogonal, not the same thing", loc="left", fontsize=10.6, fontweight="bold")
    ax.set_xticks([0, 15, 40, 70, 90, 144])

    fig.savefig(os.path.join(HERE, "2d_check.pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(HERE, "2d_check.png"), dpi=300, bbox_inches="tight")
    print("Wrote 2d_check.pdf and 2d_check.png")


if __name__ == "__main__":
    main()
