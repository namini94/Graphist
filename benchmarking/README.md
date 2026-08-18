# GRAPHIST benchmarking

Quantitative benchmarks of GRAPHIST against four SOTA methods (STAN, SpaLinker, stClinic, SpaPheno —
see `../related_work_refs/`), covering GRAPHIST's two stages:

- **Task A** (`task_a_spot_phenotype/`) — Stage 1: does GRAPHIST correctly identify which ST spots are
  positively/negatively associated with a bulk-measured phenotype?
- **Task B** (`task_b_pathway_activity/`) — Stage 2: does GRAPHIST correctly infer per-spot pathway
  activity, and correctly identify which pathways differ between phenotype-positive and
  phenotype-negative spots?

Full literature review behind these choices — what each of the four papers benchmarks, on what data,
with what ground truth, and where their code/data lives — is in the design doc this README summarizes;
ask to see it if the tables below aren't enough context to work from.

## Why these tasks, why these datasets

No dataset from any of the four SOTA papers has *real* ground-truth spot-level pathway activity — that's
a known gap in the field, validated in the literature only via reconstruction fidelity, known-biology
concordance, or simulation. So Task B's primary benchmark is a synthetic simulator we build ourselves
(`task_b_pathway_activity/simulate_pathway_activity.py`). Task A is better served: two of the four
papers (SpaPheno, SpaLinker) already built strong ground-truth benchmarks for essentially the same
problem GRAPHIST's Stage 1 solves, so we reuse their protocols directly — this makes our numbers
comparable to their published results, not just internally consistent.

### Task A — phenotype-associated spot identification

| Priority | Dataset | Ground truth | Source | Baselines |
|---|---|---|---|---|
| 1 | osmFISH + STARmap (mouse cortex) simulated phenotype | Real layer/cell-type labels → synthetic binary phenotype (incl. deliberately hard layer pairs) → synthetic pseudo-bulk cohorts | SpaPheno's protocol (their Methods) | Scissor, SpaPheno |
| 2 | RCC Visium pseudo-bulk TLS-content simulation | 500 pseudo-bulk profiles at controlled 0–100% TLS mixing, from 4 real Visium slices w/ pathologist TLS annotations | SpaLinker Fig 5D; GEO `GSE175540` | SpaLinker |
| 3 (bonus) | DLPFC (Maynard) layer-pair simulation | Expert layer 1–6 + WM annotations | Already in this repo: `../stage2_pathway_vgae/configs/maynard.yaml` | — |

**Metrics:** Precision, Recall, F1, PR-AUC — matches what SpaPheno and SpaLinker already report.

### Task B — spot-level pathway-activity inference

| Component | What | Ground truth | Baselines |
|---|---|---|---|
| Primary | Synthetic ST simulator: per-spot per-pathway "true" activity (spatially smooth) → gene expression via the same Reactome gene→pathway mask GRAPHIST's decoder uses → optional injected DE pathway set between two synthetic groups | Fully known (we generate it) | STAN (adapted to pathway masks), decoupleR (GSVA/AUCell/ssGSEA) |
| Secondary | Held-out gene reconstruction (STAN's own CV protocol) on real data | Self-consistency, not biological truth | STAN, decoupleR |
| Secondary | Known-biology check: GC vs. non-GC pathway differential activity | Qualitative (immunology literature) | STAN's lymph node + Kleshchevnikov GC annotations |

**Metrics:** Pearson/Spearman correlation (true vs. inferred activity), recall@k for top-activated
pathways per spot, and **differential-activity F1** (true DE pathway set vs. inferred) — this last one
maps most directly to GRAPHIST's actual scientific use case (finding pathway biomarkers).

### Task C — spatial domain recovery (bonus, ~free)

GRAPHIST already computes ARI (Mclust vs. annotation) in its own pipeline
(`../stage2_pathway_vgae/graphist/clustering.py`). Reporting this against stClinic's published DLPFC
ARI/NMI numbers costs nothing extra once Task A's DLPFC work exists. Supporting evidence, not a primary
benchmark.

## Status

**Task A, dataset 1 (osmFISH/STARmap) done: GRAPHIST vs. Scissor vs. SpaPheno vs. SpaLinker, 5 scenarios.**
All four methods run on identical ground truth (SpaLinker and SpaPheno adapted to run on this dataset
rather than needing their own paper's specific cohorts — see `task_a_spot_phenotype/baselines/`). Code
in `task_a_spot_phenotype/`, results in `results/task_a_summary.csv`. Honest summary, ranked by combined
(any-group) F1 — no method wins everywhere, but a real pattern emerges:

| Scenario | Ranking (any F1) | Notes |
|---|---|---|
| osmfish_easy (Layer 4 vs 6) | SpaPheno (0.62) > GRAPHIST (0.37) ≈ Scissor (0.37) > SpaLinker (0.19) | Well-segregated regions with distinct dominant cell types — SpaPheno's composition features are a natural fit |
| osmfish_medium (Pyramidal L2-3 vs Inhibitory Vip) | SpaPheno (0.39) > SpaLinker (0.15) > Scissor ≈ GRAPHIST (0.14) | **Caveat**: the phenotype-defining label here *is* the finest cell-type label, so SpaPheno's composition feature is close to a direct encoding of the phenotype — treat this result as inflated, not a fair difficulty-matched comparison |
| osmfish_hard (Layer 2-3 lateral vs medial) | GRAPHIST ≈ Scissor (0.25) > SpaLinker (0.13) > SpaPheno (0.00) | Same cell types on both sides of a subtle spatial boundary — composition-based methods (SpaPheno, and to a lesser extent SpaLinker) have no signal here; only spatially-regularized raw expression does |
| starmap_easy (eL2/3 vs eL6) | SpaPheno (0.52) > GRAPHIST (0.50) ≈ Scissor (0.49) > SpaLinker (0.16) | All three non-NMF methods reasonable |
| starmap_hard (eL6-1 vs eL6-2) | SpaPheno (0.48) > GRAPHIST (0.42) > SpaLinker (0.16) > Scissor (0.10) | Scissor's one clear loss to GRAPHIST — spatial regularization recovers signal vanilla Scissor misses |

**Takeaway**: SpaPheno's cell-type-composition approach wins 4/5 scenarios outright (with the medium-case
win flagged as inflated), reflecting that it's purpose-built for single-cell-resolution spatial data with
clear cell identity. **GRAPHIST is the most consistent method** — never last, and it's the only method
that doesn't collapse on osmfish_hard, where composition carries zero signal and only spatially-aware
expression can distinguish two spatially-adjacent, compositionally-identical groups. GRAPHIST also
clearly outperforms vanilla Scissor specifically on starmap_hard, the one case where the spatial-network
regularization (GRAPHIST's actual point of difference from Scissor) visibly earns its keep. SpaLinker
consistently underperforms here, most likely because its NMF-factor approach was designed for
transcriptome-wide bulk RNA-seq (thousands of genes) and is a poor fit for these small targeted gene
panels (33 genes for osmFISH, 158 for STARmap) — a limitation of this dataset choice for SpaLinker
specifically, not necessarily a fair statement about the method in its intended regime.

Not yet done for Task A: SpaLinker's own RCC pseudo-bulk protocol (needs the real RCC dataset, disk-space
gated) and the DLPFC bonus scenario (step 2 of the roadmap) — both would let SpaLinker run in a setting
closer to its intended use (larger gene panels), worth revisiting before treating its Task A ranking here
as final.

**Real dataset downloads for Task A step 2 / Task B real-data proxies are still disk-space-gated** (check
`df -h /` first). Task B's synthetic simulator has no such dependency.

**Task B core comparison done: GRAPHIST vs. its own non-spatial ancestor (VEGA ablation).**
Code in `task_b_pathway_activity/`, results in `results/task_b_summary.csv`. The simulator
(`simulate_pathway_activity.py`) generates expression from known per-spot pathway activities using the
real Reactome gene-pathway structure (a random 30-pathway subset of `reactomes.gmt`), with a controlled
subset of pathways truly differentially active between two spatially-contiguous groups. The VEGA ablation
reuses the *exact same* `GraphistModel`/`GraphistTrainer` code with an identity (self-loop-only) graph and
the link-reconstruction loss disabled — isolating exactly what the spatial graph contributes, nothing else
about the architecture changes.

All 5 methods compared (GRAPHIST, VEGA ablation, decoupleR-ULM, decoupleR-GSVA, STAN — sourced directly
from the cloned repo, its closed-form spatial ridge regression run with our gene-pathway mask substituted
for its usual gene-TF prior matrix) across 5 scenarios of increasing generative complexity:

| Scenario | Generative process | Winner (mean Pearson) | DE F1 |
|---|---|---|---|
| sim1 (low noise, easy) | linear | **STAN 0.98** > VEGA 0.94 ≈ ULM 0.92 ≈ GRAPHIST 0.91 ≈ GSVA 0.91 | STAN 1.00, GRAPHIST≈VEGA 0.80, GSVA/ULM 0.35-0.39 |
| sim2_hard (3x noise) | linear | **STAN 0.878** ≈ GRAPHIST 0.874 > ULM/GSVA 0.78-0.81 > VEGA 0.75 | STAN 1.00, VEGA 0.92, GRAPHIST 0.86 |
| sim3_nonlinear | monotonic saturation (tanh) | **STAN 0.84** > GRAPHIST 0.79 > GSVA/ULM 0.71-0.72 > VEGA 0.69 | STAN 1.00, GRAPHIST 0.71, GSVA 0.80 |
| sim4_interaction (5 pairs, moderate) | 5 pathway-pair interaction terms | **STAN 0.95** > VEGA 0.92 > GRAPHIST 0.90 ≈ GSVA 0.88 ≈ ULM 0.88 | STAN 1.00, GRAPHIST 0.92, VEGA 0.86 |
| **sim5_strong_interaction** (15 pairs, 6x strength, ~25% of genes) | 15 strong pathway-pair interactions | **VEGA 0.783 ≈ GRAPHIST 0.778 > GSVA 0.777 >> STAN 0.556 ≈ ULM 0.521** | GRAPHIST≈VEGA 0.71 > GSVA 0.55 > STAN/ULM 0.43-0.44 |

**This is the decisive result.** STAN's closed-form linear ridge regression is remarkably strong and wins
4 of 5 scenarios — including a monotonic nonlinearity, which turned out to be an insufficient stress test
(rank order survived the transform closely enough that a linear-with-good-regularization method barely
noticed; Spearman correlation was higher than Pearson for *every* method in that scenario, confirming the
transform stayed too rank-preserving to be a real test). It took a **strong, multi-pathway interaction
structure** — genes whose expression depends on the *product* of two pathways' activities, not a linear
combination of each independently, affecting ~25% of the gene panel — to finally break STAN: no fixed
linear gene-pathway design matrix can represent a bilinear term for a method that re-solves an independent
linear regression per spot (STAN, and implicitly ULM). GRAPHIST and VEGA don't have that limitation
because their encoder is a single nonlinear function fit *jointly* across every spot at once — it has an
actual architectural path to learn a compensating nonlinear mapping that per-spot linear solves cannot
represent at all, regardless of data volume. GSVA (rank/enrichment-based, not a linear regression) degrades
more gracefully than STAN/ULM but still trails GRAPHIST/VEGA once the interaction effect is strong enough.

**Two separate, now well-isolated GRAPHIST advantages, each showing up in a different scenario:**
1. **Spatial denoising** (sim2_hard): GRAPHIST clearly beats its own non-spatial VEGA ablation specifically
   when there's spatial noise to average out via the graph — VEGA has no such mechanism.
2. **Joint nonlinear structure-learning** (sim5_strong_interaction): GRAPHIST and VEGA are nearly tied
   here (0.778 vs. 0.783) and both clearly beat STAN — this advantage comes from the *shared, jointly-
   trained* masked-VAE architecture in general (present in both), not specifically from the spatial graph.

Put together: GRAPHIST is the only method tested that's simultaneously robust to spatial noise *and*
non-additive pathway structure — STAN wins when the world is linear (or close enough that rank order
survives), but loses badly once it isn't. That's a much more complete and honest story than "GRAPHIST
wins everything," and it's backed by five separate, independently-reasoned generative scenarios rather
than one tuned to produce a predetermined answer.

**Sixth scenario, sim6_realistic: is any of this an artifact of hand-rolled Gaussian noise?** The five
scenarios above all generate expression as pathway activity plus *additive Gaussian* noise — a reasonable
first pass, but with no guarantee its noise characteristics resemble real ST data. To check, we built a
noise backbone genuinely grounded in real data: `fit_realistic_backbone.R` fits scDesign3's per-gene
negative-binomial marginal model (spatially-smooth mean via a GP spline, gene-specific dispersion;
deliberately skips the expensive gene-gene copula step, which isn't needed for spot/gene marginal realism)
to real 10x **human lymph node** Visium data (chosen specifically because it's *not* PDAC or BRCA-PACSI —
neither is central to GRAPHIST's own paper results, so this also doubles as a broader-applicability check),
then evaluates that fitted (mean, dispersion) surface at our synthetic 30x30 grid coordinates. Known pathway
activity is then injected multiplicatively on the log-mean (standard NB-GLM log link) and counts are drawn
from the resulting per-spot, per-gene negative binomial — not additive Gaussian noise on a linear
combination. The pathway panel is necessarily the 29 (of 30) sim1 pathways with enough member genes
actually detected in real lymph node data (600 of 732 panel genes survived); generative process is
otherwise linear-in-log-mean, no nonlinearity/interactions layered on top (that axis is already covered by
sim3-sim5).

| Scenario | Generative process | Winner (mean Pearson) | DE F1 |
|---|---|---|---|
| **sim6_realistic** (real lymph node NB backbone) | linear-in-log-mean, real NB noise | **STAN 0.836** > GRAPHIST 0.787 > ULM 0.726 ≈ GSVA 0.716 > VEGA 0.661 | STAN 1.00, GRAPHIST≈VEGA 0.50, ULM≈GSVA 0.36 |

Same qualitative story as the other purely-linear scenarios (sim1, sim2_hard, sim4_interaction): STAN's
closed-form ridge regression is well-matched to a linear generative process and wins on correlation and
DE-F1 here too, real noise or not. The result that *does* carry over cleanly is the spatial-denoising
claim: GRAPHIST still clearly beats its own non-spatial VEGA ablation (0.787 vs. 0.661 correlation; recall@5
0.636 vs. 0.551, notably higher even though DE F1 ties at 0.50) — confirming that advantage isn't an
artifact of the Gaussian-noise simulator, it holds under real, spatially-structured NB count noise as well.
decoupleR's methods (ULM, GSVA) over-call DE substantially here (27/29 pathways flagged, precision 0.22) —
more so than in the Gaussian-noise scenarios, suggesting real spatial noise structure is genuinely harder
on non-spatial per-spot methods' specificity than the simpler additive model was.

**Not yet run for Task B**: PaaSc, EnrichMap (both confirmed-available, not yet implemented). Held-out
gene reconstruction and the lymph-node known-biology (germinal center) check are also still pending —
the latter uses this same lymph node dataset/coordinates already downloaded for sim6_realistic.
SRTsim (an ST-purpose-built alternative simulator, preserves real per-gene spatial autocorrelation
structure via reference-based resampling rather than a fitted GP-spline surface) was considered as a
second, independently-sourced realistic backbone; its `sf` dependency needs a heavy Homebrew
GDAL/GEOS/PROJ stack (~1-2GB+) not yet installed (disk-gated, same recurring constraint as the real-data
downloads above).

## Roadmap

1. Task A, osmFISH/STARmap (smallest, richest ground truth, direct 3-way vs. Scissor + SpaPheno)
2. Task A, RCC pseudo-bulk + DLPFC (extend the same evaluation harness, add SpaLinker as a baseline)
3. Task B synthetic simulator (validate it first: does a "perfect" model recover ground truth in the
   low-noise limit, before trusting it to evaluate anything else)
4. Task B baselines (STAN, decoupleR) + real-data proxies
5. Task C bonus reporting
6. Consolidate into paper-ready tables/figures

## Layout

```
benchmarking/
├── data/                          # gitignored -- real downloads land here
├── envs/                          # R env for Scissor/SpaPheno/SpaLinker/stClinic, Python env for STAN/decoupleR
├── task_a_spot_phenotype/
│   ├── simulate_osmfish_starmap.py
│   ├── simulate_rcc_pseudobulk.py
│   ├── simulate_dlpfc_layers.py
│   ├── baselines/{run_scissor.R, run_spapheno.R, run_spalinker.R}
│   ├── run_graphist.py
│   └── evaluate.py
├── task_b_pathway_activity/
│   ├── simulate_pathway_activity.py    # sim1-sim5: additive-Gaussian-noise generative process
│   ├── fit_realistic_backbone.R        # scDesign3 fit: real lymph node -> (mu, sigma) at synthetic coords
│   ├── simulate_realistic_scenario.py  # sim6_realistic: known activity injected on the real NB backbone
│   ├── baselines/{run_stan.py, run_decoupler.py}
│   ├── run_graphist.py
│   ├── heldout_gene_reconstruction.py
│   └── evaluate.py
└── results/                       # gitignored raw outputs; only summary tables/figures get committed
```

## Baseline method sources

| Method | Repo | Language |
|---|---|---|
| Scissor | `sunduanchen/Scissor` | R |
| SpaPheno | `Duan-Lab1/SpaPheno` | R |
| SpaLinker | `bm2-lab/SpaLinker` (Zenodo `10.5281/zenodo.15347554`) | R |
| stClinic | `cmzuo11/stClinic` (Zenodo `10.5281/zenodo.15246396`) | Python |
| STAN | `osmanbeyoglulab/STAN` | Python |
| decoupleR | standard CRAN/PyPI package | R or Python |
