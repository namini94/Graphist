"""Publication figure: does the "phenotype-differentiated tissue has more multi-pathway
crosstalk than generic tissue" finding generalize across tissues, not just the one lymph
node dataset? Summarizes estimate_real_crosstalk.py run across three independent real 10x
Visium datasets spanning three different tissue/disease contexts:
  - lymph node: germinal center (GC) vs. non-GC          [immune biology]
  - BRCA-COMMOT: Tumor vs. Healthy                        [cancer biology]
  - DLPFC/Maynard: cortical gray matter (Layers) vs. WM   [neurobiology]

For DLPFC, WM was the statistically well-powered side to test directly (background pool
must be >= phenotype-group size for sound bootstrap resampling -- oversampling a small
pool up to a much larger target duplicates points and corrupts the significance-count
metric; caught via exactly this happening on a first attempt). Sign is flipped for display
so all three tissues share one convention: positive = the biologically differentiated/
complex region (GC, Tumor, gray matter) shows MORE crosstalk than its generic counterpart.
DLPFC is shown for both the immune-themed panel (reused verbatim from the other two
tissues) and a brain-appropriate neural/synaptic panel (curated the same way
PHENOTYPE_PROGRAM was, for fairness) -- both agree.

Forest-plot style: one row per test, two z-scores (crosstalk strength, crosstalk breadth)
each with its own marker, a zero reference line, "favors hypothesis" annotated.
"""
import os

import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(__file__)

# (label, z_strength, z_breadth, sign) -- sign flips DLPFC's WM-as-tested z back to the
# shared "positive = differentiated region has more crosstalk" convention.
ROWS = [
    ("Lymph node\nGC vs. non-GC\n(immune panel)", 5.99, 3.67, 1),
    ("BRCA-COMMOT\nTumor vs. Healthy\n(immune panel)", 7.79, 2.17, 1),
    ("DLPFC\nGray matter vs. WM\n(immune panel)", -6.75, -6.82, -1),
    ("DLPFC\nGray matter vs. WM\n(neural panel)", -3.87, -4.65, -1),
]

STRENGTH_COLOR = "#2a78d6"
BREADTH_COLOR = "#eb6834"

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
    labels = [r[0] for r in ROWS]
    z_strength = np.array([r[3] * r[1] for r in ROWS])
    z_breadth = np.array([r[3] * r[2] for r in ROWS])
    y = np.arange(len(ROWS))[::-1]

    fig, ax = plt.subplots(figsize=(7.4, 3.6), constrained_layout=True)

    ax.axvspan(0, max(z_strength.max(), z_breadth.max()) + 2, color=STRENGTH_COLOR, alpha=0.045, zorder=0)
    ax.axvline(0, color="#888888", linewidth=0.9, zorder=1)

    offset = 0.16
    ax.scatter(z_strength, y + offset, s=64, color=STRENGTH_COLOR, zorder=4,
               edgecolor="white", linewidth=0.7, label="Crosstalk strength (mean |r|)")
    ax.scatter(z_breadth, y - offset, s=64, color=BREADTH_COLOR, zorder=4, marker="D",
               edgecolor="white", linewidth=0.7, label="Crosstalk breadth (# significant pairs)")
    for yi, zs, zb in zip(y, z_strength, z_breadth):
        ax.plot([min(zs, zb) - 0, max(zs, zb)], [yi, yi], color="#cccccc", linewidth=6, zorder=2, solid_capstyle="round")

    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=9.2)
    ax.set_xlabel("z-score vs. null (random same-sized generic-tissue subsets)\n"
                  "positive = the phenotype-differentiated region shows more crosstalk")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.tick_params(length=0)
    ax.set_xlim(min(z_strength.min(), z_breadth.min()) - 1.5, max(z_strength.max(), z_breadth.max()) + 1.5)

    # The shaded region (positive = favors hypothesis) plus the explicit axis-label wording
    # already convey directionality clearly; a redundant "favors/against" text overlay was
    # cramped against the data at this figure size, so it's dropped rather than forced in.
    ax.legend(loc="upper left", bbox_to_anchor=(0.0, 1.0), frameon=False, fontsize=8.6, handletextpad=0.6)
    fig.suptitle("Phenotype-differentiated tissue shows more multi-pathway crosstalk\nacross three independent tissues",
                 x=0.01, ha="left", fontsize=11, fontweight="bold", y=1.12)

    fig.savefig(os.path.join(HERE, "multi_tissue_crosstalk.pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(HERE, "multi_tissue_crosstalk.png"), dpi=300, bbox_inches="tight")
    print("Wrote multi_tissue_crosstalk.pdf and multi_tissue_crosstalk.png")


if __name__ == "__main__":
    main()
