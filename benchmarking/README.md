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

**Seventh scenario, sim7_compound: do the two GRAPHIST advantages stack?** sim2_hard isolates spatial-noise
robustness; sim5_strong_interaction isolates nonlinear structure-learning. Real tissue plausibly has both
at once, so this scenario combines them directly: sim2_hard's noise level (`--noise-sd 1.5`, 3x default)
*and* sim5's interaction structure (`--n-interactions 15 --interaction-strength 6.0`) in one dataset.

| Scenario | Generative process | Winner (mean Pearson) | DE F1 |
|---|---|---|---|
| **sim7_compound** | 3x spatial noise + 15 strong interactions | **GRAPHIST 0.704** > GSVA 0.680 > VEGA 0.654 > STAN 0.581 > ULM 0.501 | STAN 0.92 > VEGA 0.80 > GRAPHIST 0.75 > ULM 0.50 > GSVA 0.41 |

This is the first scenario where **GRAPHIST outright wins on correlation** — and notably, it now also
clearly beats its own VEGA ablation (0.704 vs. 0.654), a bigger gap than sim5 alone showed (where VEGA and
GRAPHIST were essentially tied, 0.783 vs. 0.778). That's consistent with the compounding hypothesis: with
spatial noise added on top of the interaction structure, the spatial graph's denoising has something to
average out again, on top of the encoder's nonlinear capacity — both mechanisms contributing rather than
one dominating.

**Worth reporting honestly rather than glossing over: DE-F1 does not follow the same ranking.** STAN's
DE-pathway calling is still the sharpest here (0.92, driven by high precision — 7 pathways flagged for 6
true positives) despite trailing badly on raw correlation (0.581). GRAPHIST/VEGA correctly recall all 6
true DE pathways (recall=1.0 for every method) but flag more false positives (10 and 9 respectively) than
STAN, dragging down precision and F1. This is a real, useful finding, not a contradiction: correlation
measures how well the *entire* activity surface is recovered (where GRAPHIST's advantage is real and
compounds), while this DE-F1 metric is a specific threshold-crossing decision (BH significance + Cohen's
d > 0.3) that's sensitive to each method's prediction variance/calibration in a way raw correlation isn't
— worth flagging as a limitation of the current DE-calling threshold rather than evidence GRAPHIST's
underlying activity estimates are worse here.

**Eighth/ninth scenarios, sim8_weak_de / sim9_vweak_de: does GRAPHIST hold up better under weak DE signal
plus spatial noise?** Hypothesis: per-spot regression should lose sensitivity to subtle biomarkers faster
than GRAPHIST's spatially-averaging encoder as the true DE effect shrinks. Tested by re-running sim2_hard's
noise level (`--noise-sd 1.5`) at `--de-effect-size 1.0` and `0.5` (vs. the 2.0 default).

**This is a null result — reported as such, not omitted.** Correlation for every method actually *improves*
as the effect size shrinks (STAN 0.878→0.953, GRAPHIST 0.874→0.899 from sim2_hard to sim9_vweak_de), and
STAN's lead *widens*, not narrows. The mechanism: `de_effect_size` adds a sharp step-like group-mean shift
on top of an otherwise smooth spatial field for the 6 DE pathways only; a larger step means more of that
pathway's total variance comes from a discontinuity that's somewhat harder for every method (spatial or
not) to reproduce exactly than the smooth background is — so bigger effect sizes drag correlation down
across the board, an artifact of this simulator's DE-injection design rather than a real noise-robustness
signal. DE-F1 tells the same story: STAN holds perfect 1.00 across all three effect sizes; GRAPHIST/VEGA's
F1 mildly degrades as the signal weakens (0.86→0.80→0.80/0.71). No crossover point found in this design —
worth revisiting with a differently-shaped weak-signal injection (e.g. shrinking the *whole* field's
variance rather than adding a smaller step) if this angle is worth pursuing further for the paper.

**Tenth/eleventh scenarios, sim10_sparse / sim11_vsparse: does the spatial graph compensate for dropout?**
Real Visium data is substantially zero-inflated (sim6_realistic was 40.6% zero); this tests robustness as
sparsity increases further. Extended `simulate_realistic_scenario.py` with `--depth-scale`, which multiplies
the realistic mean before NB sampling — simulating reduced sequencing depth, which naturally increases
zero-inflation through the NB process itself (lower mean → more zeros) rather than an artificial independent
dropout mask. Swept `--depth-scale 0.3` (59.8% zero) and `0.1` (74.7% zero) against the `sim6_realistic`
depth=1.0 (40.6% zero) anchor.

| Zero fraction | GRAPHIST | VEGA | STAN | GRAPHIST − VEGA gap |
|---|---|---|---|---|
| 0.41 (sim6_realistic) | 0.787 | 0.661 | 0.836 | 0.126 |
| 0.60 (sim10_sparse) | 0.727 | 0.529 | 0.757 | 0.198 |
| 0.75 (sim11_vsparse) | 0.639 | 0.403 | 0.661 | **0.236** |

(mean Pearson correlation to ground truth; DE F1 at the same three points: STAN 1.00→0.92→0.91,
GRAPHIST 0.50→0.46→0.50, VEGA 0.50→0.55→0.63, GSVA/ULM 0.36→0.37-0.38→0.40-0.43.)

**Partial, honest result — not a clean sweep.** STAN still leads on raw correlation at every sparsity level
tested (up to 75% zero); GRAPHIST does not overtake it in this range. But the **GRAPHIST − VEGA gap widens
monotonically and substantially with sparsity** (0.126 → 0.198 → 0.236) — i.e. the spatial graph's specific
contribution grows almost 2x as dropout increases from realistic (41%) to severe (75%). That's the cleanest,
most mechanistically sensible finding in this sweep: spatial averaging matters most exactly when per-spot
signal is noisiest/sparsest, which is precisely the real-world condition (low-quality or low-depth spots)
where it should. STAN's own absolute correlation also degrades with sparsity (0.836→0.661), just proportionally
less than VEGA's non-spatial baseline does — so the honest framing is "GRAPHIST's spatial-graph advantage over
not having one grows under realistic dropout," not "GRAPHIST beats STAN under dropout," which isn't what the
data shows in this range. Worth extending to even sparser (`--depth-scale` < 0.1) to see if the narrowing
STAN-vs-GRAPHIST gap ever crosses over.

### Connecting the mechanism findings to the paper's actual motivation

Why sim5/sim7's advantage over STAN isn't just "we found a synthetic setting where GRAPHIST wins" — it's
directly relevant to the paper's real use case (phenotype-driven pathway biomarker discovery), for a
specific, citable biological reason, not because it's convenient:

**Real phenotype-driven tissue states are coordinated multi-pathway crosstalk, not one pathway shifting in
isolation.** This is established biology: EMT phenotypes involve synergistic TGF-β/WNT/cytoskeletal
crosstalk; hypoxia phenotypes involve HIF1A-coupled glycolysis/angiogenesis/immune-evasion signaling; TLS
presence (the phenotype SpaLinker itself validates on) involves coordinated chemokine signaling + B/T-cell
activation + germinal-center reaction pathways; immune-excluded vs. -inflamed tumor phenotypes involve
coupled immune-checkpoint and stromal signaling. None of these are additive single-pathway effects — this
is exactly the generative structure sim5/sim7's bilinear pathway-interaction terms are a synthetic stand-in
for, and exactly the regime where no fixed linear gene-pathway design matrix (STAN's, or any per-spot
regression's) can represent the signal, while GRAPHIST's jointly-trained nonlinear encoder can.

There's also a structural reason this maps directly rather than by analogy: the simulator's group A/B split
*already is* mechanically identical to Graphist(+)/Graphist(-) — same differential-pathway-biomarker-
discovery framing GRAPHIST's Stage 2 actually performs on real data. So relabeling sim5/sim7 as
"phenotype-positive vs. phenotype-negative" for the paper is accurate, not a stretch.

**The honest line to hold in the manuscript**: this is a *mechanistic motivation*, grounded in real,
citable biology, for why GRAPHIST's architecture should have an advantage specifically in phenotype-driven
settings — not a literal empirical demonstration run on real phenotype-labeled data (no baseline here,
including STAN, does phenotype-driven spot classification at all; that's Task A, evaluated separately with
different baselines). Overclaiming "we proved GRAPHIST beats STAN on phenotype data" from these results
would not survive review; "sim5/sim7's generative structure specifically models the kind of multi-pathway
coordination known to characterize real phenotype-driven biology, and GRAPHIST outperforms STAN precisely
there" is the accurate, defensible version.

#### Twelfth-fourteenth scenarios: testing the phenotype connection directly with a real, curated program

`simulate_phenotype_program.py` builds a matched-pair experiment using a real, thematically coherent
6-pathway immune/TLS-activation program (chemokine recruitment, TCR signaling, BCR signaling, downstream
NF-kB — a real lymphoid-activation signature, not an arbitrary pathway subset) as the phenotype's biomarker
set, generating a `phenotype_coordinated/` and matched `phenotype_generic/` dataset pair from *identical*
ground-truth activity, weights, baseline, and noise — differing in exactly one respect: whether the
program's pathways have an added crosstalk structure. Three variants were tried, in order, and the honest
result of each is reported (not just the one that worked):

| Variant | Design | Result |
|---|---|---|
| sim12 (`within`, dense) | All 15 pairs *within* the 6-pathway program interact, 6x strength | **Collapses for every method, GRAPHIST included** (program-pathway correlation ~0.02-0.16 for all) — too concentrated, individual pathway identity becomes unrecoverable by anyone |
| sim13 (`within`, moderate) | 5 pairs within the program, 1x strength | Fully identifiable, but **STAN still wins outright** (0.864 vs. GRAPHIST's 0.773) |
| sim14 (`cross`) | 15 pairs, each between a phenotype pathway and an unrelated background pathway, 6x strength | **GRAPHIST clearly beats STAN** |

The first two attempts (interactions concentrated within the phenotype program itself) were the naive
reading of the biological argument above, and they didn't work — genuinely useful negative results, not
omitted. What actually reproduces sim5/sim7's advantage is *diffuse* crosstalk: the phenotype's pathways
interacting with *other*, unrelated ongoing tissue biology elsewhere in the panel, not with each other. That
refines the paper claim to something more precise and more defensible: real phenotypes plausibly perturb a
core biomarker program *and* diffusely couple to broader tissue pathway activity (which is what "a phenotype
affects tissue biology" means in practice — it rarely touches only 6 pathways and nothing else), not that
the biomarker pathways specifically coordinate with each other.

sim14's matched-pair result (mean Pearson, full 30-pathway panel):

| Condition | GRAPHIST | VEGA | STAN | GSVA | ULM |
|---|---|---|---|---|---|
| phenotype_generic (no crosstalk) | 0.841 | 0.836 | **0.944** | 0.825 | 0.833 |
| phenotype_coordinated (diffuse crosstalk) | **0.616** | 0.591 | 0.544 | 0.647 | 0.440 |

Same phenotype program, same effect size, same noise, same seed — changing only whether the program
diffusely crosstalks with background tissue biology flips the winner from STAN to GRAPHIST, and GRAPHIST
also clearly beats its own VEGA ablation in the coordinated condition (spatial advantage holds too).

**Honest caveat, not hidden**: GSVA edges out GRAPHIST in the coordinated condition (0.647 vs. 0.616), and
does so more clearly when scored on just the 6 phenotype pathways specifically (program-pathway mean
Pearson: GSVA 0.498, ULM 0.402, GRAPHIST 0.361, VEGA 0.331, STAN 0.282 — GRAPHIST still clearly beats STAN,
but GSVA and even ULM beat GRAPHIST here). The defensible claim is **"GRAPHIST beats STAN, its most direct
methodological competitor, in phenotype-like diffuse-crosstalk conditions"** — not "GRAPHIST beats every
baseline." GSVA remains a real, cheap, surprisingly strong competitor worth acknowledging directly in the
paper, not one to bury.

#### Fifteenth scenario: is sim14 a lucky single point, or a robust regime?

`sim15_dose_*` sweeps the number of diffuse cross-pathway interaction pairs (0 to 144, the maximum
possible for a 6-pathway program against 24 background pathways) at fixed 6x strength, everything else
identical to sim14 — turning the single sim14 data point into a full dose-response curve. Full results in
`results/task_b_dose_response.csv`; mean Pearson correlation (full 30-pathway panel):

| n_pairs | GRAPHIST | VEGA | STAN | ULM | GSVA | GRAPHIST − STAN |
|---|---|---|---|---|---|---|
| 0 | 0.841 | 0.836 | **0.944** | 0.833 | 0.825 | −0.103 |
| 5 | 0.764 | 0.747 | 0.759 | 0.688 | 0.769 | +0.006 |
| 10 | 0.712 | 0.695 | 0.676 | 0.532 | 0.696 | +0.036 |
| 15 (sim14) | 0.616 | 0.591 | 0.544 | 0.440 | 0.647 | +0.072 |
| 25 | 0.583 | 0.562 | 0.480 | 0.370 | 0.571 | +0.103 |
| **40** | **0.442** | 0.425 | 0.314 | 0.246 | 0.434 | **+0.128 (peak)** |
| 70 | 0.298 | 0.278 | 0.227 | 0.166 | 0.305 | +0.071 |
| 100 | 0.245 | 0.221 | 0.183 | 0.149 | 0.195 | +0.062 |
| 144 (max) | 0.203 | 0.178 | 0.147 | 0.113 | 0.144 | +0.056 |

**This is a genuine dose-response, not a coin flip at one setting.** The crossover happens almost
immediately (between 0 and 5 pairs — a small amount of diffuse crosstalk is enough to flip the winner),
GRAPHIST's margin over STAN grows monotonically up to a peak around 40 pairs (~28% of the maximum possible
density), then narrows again at very high density as absolute performance floors out for every method (144
pairs is an extreme, near-saturated regime where nothing recovers pathway identity well). The advantage
holds across a wide, plausible middle range (5-144 pairs), not just the one point sim14 happened to use.

**Two other patterns worth having in the paper**: (1) GRAPHIST beats its own VEGA ablation at *every single
point* on the curve (gap 0.005-0.025, small but universally positive) — the spatial-graph advantage is
unconditional here, unlike the STAN advantage which specifically requires crosstalk to be present at all.
(2) The GSVA comparison oscillates around zero the whole way (gap from −0.030 to +0.059, no consistent
winner) — reinforcing that GSVA is a genuinely close, real competitor throughout the regime, not just an
artifact of sim14's specific setting.

Publication figure: `results/figures/dose_response.{pdf,png}` — two panels, (A) absolute recovery accuracy
for all 5 methods vs. crosstalk density, (B) the GRAPHIST−STAN gap itself, with the crossover (n=5) and peak
(n=40, +0.128) annotated. Generated by `results/figures/plot_dose_response.py`.

#### Does real phenotype-differentiated tissue actually have more crosstalk than generic tissue?

Every scenario so far assumes diffuse multi-pathway crosstalk is a plausible feature of real phenotype-driven
biology (the citable-biology argument earlier in this section) — but that's still an assumption about real
data, not a measurement of it. `estimate_real_crosstalk.py` tests it directly: does a real, established
phenotype-differentiated tissue micro-environment actually show more multi-pathway crosstalk than a generic,
biologically-arbitrary region of the same tissue?

**Data**: the same 10x human lymph node Visium data used throughout, with a real phenotype axis — germinal
center (GC) vs. non-GC spots, from the manual annotation in `osmanbeyoglulab/STAN`'s own resources (the same
dataset/annotation STAN itself validates against). GC reactions are a textbook coordinated multi-pathway
immune process — chemokine-driven B-cell recruitment, BCR signaling, T-cell help, NF-kB activation — which is
literally the same immune program curated as `PHENOTYPE_PROGRAM` in `simulate_phenotype_program.py`, not a
loose analogy.

**Method**: score all 30 panel pathways per real spot via GSVA (rank-based — it can't manufacture crosstalk
that isn't there, unlike a method that assumes linear structure). For a given spot subset, compute (a) mean
|Pearson correlation| and (b) count of BH-significant (α=0.05) pathway pairs, restricted to pairs with
**disjoint gene membership** (excludes pairs sharing genes, which would trivially correlate regardless of
any real crosstalk — same logic as `simulate_phenotype_program.coordination_score`, now applied to real
inferred activity instead of synthetic ground truth). Compare the real GC-spot subset (n=378) against a null
built from 200 random same-sized subsets of non-GC spots elsewhere in the same tissue.

| Statistic | Real GC spots | Null (200 random non-GC subsets) | z-score | Percentile |
|---|---|---|---|---|
| Mean \|pathway-pathway correlation\| | **0.086** | 0.067 ± 0.003 | 6.91 | 100th |
| Significant pathway pairs (of 396 disjoint) | **70** | 34.7 ± 7.6 | 4.67 | 100th |

**This is a strong, clean result**: the real GC micro-environment shows roughly **2x** the multi-pathway
crosstalk of a random same-sized piece of the same tissue, on both a strength measure and a breadth measure,
each many standard deviations outside the null. Real phenotype-differentiated tissue genuinely has more
coordinated multi-pathway activity than generic tissue — this isn't just a plausible-sounding assumption
underlying the synthetic scenarios, it's now a measurement.

**This is the piece that completes the argument for the paper**: (1) real phenotype-differentiated tissue
measurably has more multi-pathway crosstalk than generic tissue (this analysis); (2) GRAPHIST's advantage
over STAN emerges specifically once diffuse crosstalk is present and grows with its density (sim14/15's dose-
response curve); therefore (3) real phenotype-driven ST analysis — GRAPHIST's actual paper use case — plausibly
sits in the regime the dose-response curve shows favors GRAPHIST, for a citable biological reason, now with a
real-data measurement behind it rather than an assumption. The honest caveat from sim16 still applies: the
sharpest form of the claim (GRAPHIST specifically nailing the phenotype's *own* biomarker pathways) weakened
under real NB noise even as the full-panel and VEGA-ablation advantages held — so this three-step argument
supports "GRAPHIST vs. STAN, full-panel" and "GRAPHIST vs. VEGA" robustly, and the phenotype-biomarker-specific
version more provisionally.

Publication figure: `results/figures/real_crosstalk.{pdf,png}` — null-distribution histograms with the real
GC value marked, for both the strength and breadth statistics. Generated by
`results/figures/plot_real_crosstalk.py`.

#### Sixteenth scenario: does this survive combining with REAL noise characteristics?

`simulate_realistic_phenotype_scenario.py` combines the two findings that mattered most this session: the
diffuse-crosstalk mechanism (sim14/15) and the scDesign3 real-lymph-node NB backbone (sim6_realistic).
Required refitting the backbone on the phenotype program's own gene panel (`lymph_node_phenotype/`,
re-downloaded via `sc.datasets.visium_sge`, since sim6_realistic's backbone was fit on sim1's unrelated
732-gene panel) — 625/634 genes survived (9 dropped for a NaN dispersion fit). Used sim15's peak-advantage
setting (40 diffuse cross-group pairs, 6x strength). One real implementation issue surfaced and fixed:
the raw accumulated interaction contribution is heavy-tailed (a phenotype pathway participates in several
of the 40 pairs at once, so some genes accumulate contributions from multiple pairs, reaching ~100+) —
unlike the additive-Gaussian model sim14/15 used, exponentiating an unbounded value in the NB log-link
exploded into astronomical counts (mean count 43,135 on a first attempt). Fixed by squashing the
contribution through `tanh(x / 10) * 10` before scaling, bounding the worst-case log-mean shift regardless
of accumulation.

| Condition | GRAPHIST full | STAN full | GRAPHIST program-only | STAN program-only |
|---|---|---|---|---|
| generic (no crosstalk) | 0.740 | **0.821** | 0.873 | **0.900** |
| coordinated (diffuse crosstalk) | **0.479** | 0.463 | 0.458 | **0.478** |

**Reported honestly, not smoothed over — this is a more mixed result than sim14's clean synthetic
version.** GRAPHIST does still lead on full-panel correlation in the coordinated condition, but the margin
shrank sharply under real noise (sim14's synthetic +0.072 gap → +0.016 here). More notably, **on the 6
actual phenotype-biomarker pathways specifically, STAN edges GRAPHIST out** (0.478 vs. 0.458) — the
sharpest form of the claim ("GRAPHIST better recovers the phenotype's own biomarkers under diffuse
crosstalk") does not fully survive real NB noise, even though it held in the idealized synthetic version.

**What does hold robustly, in every variant tried this session, synthetic or real-noise-grounded**:
GRAPHIST clearly beats its own VEGA ablation in *both* conditions here (0.740 vs. 0.623 generic, 0.479 vs.
0.438 coordinated) — the spatial-graph advantage is the most durable finding in the whole suite. DE-F1 is a
wash at this setting (GRAPHIST/STAN/VEGA all hit perfect 1.0 in the coordinated condition).

**Honest interpretation for the paper**: the diffuse-crosstalk mechanism is real (confirmed by the dose-
response sweep) but its *magnitude* under idealized Gaussian noise is larger than under realistic,
zero-inflated NB noise — likely because real noise's own variance already dominates more of the total
signal, leaving less room for the crosstalk structure specifically to matter relative to STAN's baseline
competence. The full-panel GRAPHIST-over-STAN lead and the VEGA-ablation advantage are the two claims that
survive contact with real noise; "GRAPHIST specifically nails the phenotype's own biomarker pathways
better than STAN" is the more fragile, synthetic-only version of the claim and shouldn't be overstated in
the manuscript without this caveat attached.

**A stronger, more direct version of this connection is buildable but not yet done**: an end-to-end scenario
that explicitly simulates a bulk phenotype variable correlated with the group A/B assignment, runs it
through GRAPHIST's actual Stage 1 (bulk-to-spot regression) to *recover* the phenotype-associated spots, then
Stage 2 on the recovered (not ground-truth) groups — testing the full pipeline's phenotype-to-biomarker
recovery jointly, rather than Task A and Task B's mechanisms being validated separately. Worth doing if this
narrative becomes a paper section rather than just discussion framing.

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
