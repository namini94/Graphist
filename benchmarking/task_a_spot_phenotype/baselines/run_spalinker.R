#!/usr/bin/env Rscript
# SpaLinker baseline: bulk expression + continuous phenotype -> NMF factor discovery
# -> project the phenotype-correlated factor onto ST via NMFpredict -> threshold.
#
# Ported from SpaLinker's pheno_program_NMF.R / FeatureSelection.R, sourced directly
# from the cloned repo (only the NMF fitting/prediction and correlation-based feature
# selection functions -- not the Seurat-object-coupled plotting wrapper PredNMFinST,
# which we bypass in favor of calling NMFpredict directly on a plain matrix).
#
# Pipeline (mirrors SpaLinker's Fig 5 RCC pseudo-bulk protocol, generalized to our
# osmFISH/STARmap scenarios instead of requiring the RCC dataset):
#   1. PhenoAssoFeatures(method="cor"): correlate each gene with the bulk phenotype,
#      keep significant genes as the NMF input (their own gene-selection step).
#   2. RunNMFtest: NMF-factorize the selected bulk genes x samples matrix.
#   3. Correlate each factor's H row (bulk) against the phenotype; the single most
#      correlated factor is "the phenotype factor" (sign tells us its direction).
#   4. NMFpredict: project the fixed W (gene loadings) onto ST expression to get a new
#      H matrix (factors x ST cells).
#   5. Sign-adjusted phenotype-factor row = per-cell score, thresholded with the same
#      top/bottom cutoff-fraction convention as the other three methods.
#
# Usage: Rscript run_spalinker.R <data_dir> <out_csv>

args <- commandArgs(trailingOnly = TRUE)
data_dir <- args[1]
out_path <- args[2]

spalinker_dir <- file.path(data_dir, "..", "..", "repos", "SpaLinker", "R")
suppressMessages({
  library(NMF)
  library(foreach)
  library(progressr)
})
source(file.path(spalinker_dir, "FeatureSelection.R"))
source(file.path(spalinker_dir, "pheno_program_NMF.R"))

bulk <- as.matrix(read.csv(file.path(data_dir, "bulk_expression.csv"), row.names = 1, check.names = FALSE))
phenotype_df <- read.csv(file.path(data_dir, "bulk_phenotype.csv"), row.names = 1)
phenotype <- as.numeric(phenotype_df[, 1])
names(phenotype) <- rownames(phenotype_df)
phenotype <- phenotype[colnames(bulk)]

st <- as.matrix(read.csv(file.path(data_dir, "st_expression.csv"), row.names = 1, check.names = FALSE))
common_genes <- intersect(rownames(bulk), rownames(st))
bulk <- bulk[common_genes, , drop = FALSE]
st <- st[common_genes, , drop = FALSE]

message("Selecting phenotype-correlated genes (PhenoAssoFeatures, method='cor')...")
asso <- PhenoAssoFeatures(data = bulk, phenotype = phenotype, method = "cor", verbose = FALSE)
selected_genes <- rownames(asso)[!is.na(asso$p.value) & asso$p.value < 0.1]
if (length(selected_genes) < 10) {
  message(sprintf("Only %d genes passed p<0.1; relaxing to top 50%% by |corr| instead.", length(selected_genes)))
  selected_genes <- rownames(asso)[order(-abs(asso$corr))][1:max(10, floor(nrow(asso) / 2))]
}
message(sprintf("Selected %d/%d genes for NMF.", length(selected_genes), nrow(bulk)))

bulk_nmf_input <- NMF_bulk_input(bulk[selected_genes, , drop = FALSE], dolog = FALSE)

message("Running NMF (rank=3, nrun=10)...")
nmf_fit <- RunNMFtest(expr = bulk_nmf_input, rank = 3, nrun = 10, verbose = FALSE)

H_train <- nmf_fit@fit@H
factor_corr <- apply(H_train, 1, function(h) cor(h, phenotype[colnames(H_train)]))
best_factor <- names(which.max(abs(factor_corr)))
factor_sign <- sign(factor_corr[best_factor])
message(sprintf("Phenotype factor: %s (train correlation = %.3f)", best_factor, factor_corr[best_factor]))

message("Projecting onto ST data (NMFpredict)...")
W <- nmf_fit@fit@W
st_nmf <- NMFpredict(W, st)
H_st <- st_nmf@fit@H
# NMFpredict() itself doesn't carry over factor names (that renaming normally happens in
# the Seurat-object wrapper PredNMFinST, which we bypass) -- restore them here so
# `best_factor` (identified from the training H) indexes the correct row.
rownames(H_st) <- rownames(H_train)
scores <- factor_sign * H_st[best_factor, ]
names(scores) <- colnames(H_st)
# NMFpredict silently drops all-zero-expression cells; anything missing gets a neutral score.
missing <- setdiff(colnames(st), names(scores))
if (length(missing) > 0) {
  fill <- rep(0, length(missing))
  names(fill) <- missing
  scores <- c(scores, fill)
}
scores <- scores[colnames(st)]

cutoff <- 0.2
n_cells <- length(scores)
n_side <- floor(n_cells * cutoff / 2)
ranked <- order(scores, decreasing = TRUE)
predicted <- rep(0, n_cells)
names(predicted) <- names(scores)
predicted[names(scores)[ranked[1:n_side]]] <- 1
predicted[names(scores)[rev(ranked)[1:n_side]]] <- -1

message(sprintf("Selected %d positive, %d negative (top/bottom %.0f%% by score, cutoff=%.2f)",
                 n_side, n_side, 100 * cutoff / 2, cutoff))

write.csv(data.frame(cell_id = names(predicted), predicted = predicted, score = scores, used_cutoff = cutoff),
          out_path, row.names = FALSE)
message(sprintf("Wrote predictions for %d cells to %s", length(predicted), out_path))
