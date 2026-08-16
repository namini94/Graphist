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

Design phase complete (this README + directory skeleton). Nothing has been run yet — see the phased
roadmap below. **Real dataset downloads are blocked on disk space** (check `df -h /` before starting
Task A or Task B's real-data proxies; the synthetic simulator in Task B has no such dependency and can
be built first if disk space isn't available yet).

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
│   ├── simulate_pathway_activity.py
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
