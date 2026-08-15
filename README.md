# Graphist

Code and materials for **GRAPHIST** (GRAph learning method for PHenotype-based Interpretable pathway
activity identification of Spatial Transcriptomics) — a two-stage framework that (1) links spatial
transcriptomics (ST) spots to a clinical/experimental phenotype measured on bulk RNA-seq samples, then
(2) uses a spatially-aware, pathway-masked variational graph autoencoder to quantify which gene pathways
are active in the phenotype-associated tissue regions. See [`paper/GRAPHIST_RECOMB2025.pdf`](paper/GRAPHIST_RECOMB2025.pdf)
for the full method description.

Applied to breast cancer (BRCA-PACSI, BRCA-COMMOT), pancreatic cancer (PDAC), and DLPFC (Maynard)
spatial datasets.

## Status

This repo is being reorganized out of several previously scattered project directories. The core Stage 2
model code (`stage2_pathway_vgae/`) is in the middle of being refactored from four large, near-duplicate,
hardcoded-path scripts into a clean, installable `graphist` Python package (config-driven, one shared
model implementation, tests). Until that refactor lands, `stage2_pathway_vgae/` still contains the
original per-dataset scripts as-is.

Raw datasets and large model-output artifacts are **not** included in this repo (they live locally under
`~/Documents/BulkToST/Dataset` and `~/Documents/BulkToST/Res-*`) — this repo holds code, configs, docs,
and figures only.

## Layout

- **`stage1_bulk_regression_R/`** — Stage 1: a Scissor-adapted sparse regression (R) that scores each ST
  spot for association with a bulk-sample phenotype (survival, gene knockdown, etc.), split by dataset.
- **`stage2_pathway_vgae/`** — Stage 2: the pathway-activity variational graph autoencoder (Python/PyTorch
  + PyG). Trains a GCN encoder + pathway-masked decoder per dataset; `preprocessing/` holds the scanpy
  preprocessing scripts each pipeline depends on.
- **`downstream_analysis/`** — visualization and interpretation scripts run after Stage 1/2 (spatial plots
  of Graphist+/- selections, pathway lollipop/volcano/boxplot figures, histology-morphology clustering,
  TCGA bulk data acquisition).
- **`paper/`** — the paper draft, the working figure deck, and each finalized figure (Fig2, Fig4, SFig1).
- **`related_work_refs/`** — reference PDFs for methods cited in Related Work (SpaLinker, stClinic, STAN,
  SpaPheno).
- **`generated_figures/`** — exploratory plot outputs already produced by the downstream analysis scripts.

## Requirements

No environment is pinned yet (tracked as part of the Stage 2 refactor). The Stage 2 model currently runs
against: Python 3.10, PyTorch 1.13.1, PyTorch Geometric 2.5.3, scanpy 1.9.1, anndata 0.8.0,
scikit-learn 1.2.1. Stage 1 is plain R (Scissor + standard Bioconductor/CRAN packages).
