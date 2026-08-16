#!/usr/bin/env Rscript
# Vanilla Scissor baseline: bulk expression + phenotype -> Scissor(+/-) spot calls,
# using an EXPRESSION-similarity network (Seurat SNN graph) -- Scissor's own internal
# default -- instead of GRAPHIST's spatial k-NN network. Structurally identical to
# ../run_graphist.R (same APML1 call, same penalty="Net", same alpha search) so the
# only difference between the two is the source of the network regularizer, isolating
# exactly what "spatial-awareness" buys GRAPHIST over vanilla Scissor.
#
# Note: Scissor's own high-level Scissor() wrapper is not used directly here because it
# validates `family="gaussian"` phenotypes as a small number of discrete groups (via a
# `tag` argument matching `table(phenotype)`), which doesn't fit our continuous
# pseudo-bulk mixing-fraction phenotype -- APML1 (which the wrapper itself calls
# internally) has no such restriction.
# Usage: Rscript run_scissor.R <data_dir> <out_csv> [family]

args <- commandArgs(trailingOnly = TRUE)
data_dir <- args[1]
out_path <- args[2]
family <- if (length(args) >= 3) args[3] else "gaussian"

suppressMessages(library(Scissor))
suppressMessages(library(Seurat))
suppressMessages(library(preprocessCore))
options(Seurat.object.assay.version = "v3")

bulk <- as.matrix(read.csv(file.path(data_dir, "bulk_expression.csv"), row.names = 1, check.names = FALSE))
phenotype_df <- read.csv(file.path(data_dir, "bulk_phenotype.csv"), row.names = 1)
phenotype <- as.numeric(phenotype_df[, 1])
names(phenotype) <- rownames(phenotype_df)
phenotype <- phenotype[colnames(bulk)]

st <- as.matrix(read.csv(file.path(data_dir, "st_expression.csv"), row.names = 1, check.names = FALSE))
storage.mode(st) <- "double"

common_genes <- intersect(rownames(bulk), rownames(st))
bulk <- bulk[common_genes, , drop = FALSE]
st <- st[common_genes, , drop = FALSE]

sc_obj <- CreateSeuratObject(counts = st)
sc_obj <- FindVariableFeatures(sc_obj, selection.method = "vst", verbose = FALSE)
sc_obj <- ScaleData(sc_obj, verbose = FALSE)
n_pcs <- min(10, ncol(sc_obj) - 1, nrow(sc_obj) - 1)
sc_obj <- RunPCA(sc_obj, features = VariableFeatures(sc_obj), npcs = n_pcs, verbose = FALSE)
sc_obj <- FindNeighbors(sc_obj, dims = 1:n_pcs, verbose = FALSE)
Omega <- as.matrix(sc_obj@graphs$RNA_snn)
diag(Omega) <- 0
Omega[Omega != 0] <- 1
Omega <- Omega[colnames(st), colnames(st)]

# Quantile-normalize bulk and ST together (matches Scissor's own internal preprocessing).
combined <- cbind(bulk, st)
combined_norm <- normalize.quantiles(as.matrix(combined))
rownames(combined_norm) <- rownames(combined)
colnames(combined_norm) <- colnames(combined)
bulk_norm <- combined_norm[, 1:ncol(bulk), drop = FALSE]
st_norm <- combined_norm[, (ncol(bulk) + 1):ncol(combined_norm), drop = FALSE]

X <- cor(bulk_norm, st_norm)
Y <- as.matrix(phenotype)

alpha_seq <- c(0.03, 0.1, 0.3, 0.5)
# Try cutoff=0.2 (the convention used throughout stage1_bulk_regression_R/) first; only
# relax it if nothing in alpha_seq converges under it, exactly mirroring ../run_graphist.R
# so both methods are held to the same convergence standard.
cutoff_seq <- c(0.2, 0.35, 0.5)
predicted <- rep(0, ncol(X))
scores <- rep(0, ncol(X))
names(predicted) <- colnames(X)
names(scores) <- colnames(X)
converged <- FALSE
used_cutoff <- NA
for (cutoff in cutoff_seq) {
  for (a in alpha_seq) {
    set.seed(123)
    fit0 <- tryCatch(
      APML1(X, Y, family = family, penalty = "Net", alpha = a, Omega = Omega,
            nlambda = 100, nfolds = min(10, nrow(X))),
      error = function(e) { message(sprintf("cutoff=%s alpha=%s (cv) failed: %s", cutoff, a, conditionMessage(e))); NULL }
    )
    if (is.null(fit0)) next
    fit1 <- APML1(X, Y, family = family, penalty = "Net", alpha = a, Omega = Omega, lambda = fit0$lambda.min)
    Coefs <- as.numeric(fit1$Beta)
    pos <- colnames(X)[which(Coefs > 0)]
    neg <- colnames(X)[which(Coefs < 0)]
    pct <- (length(pos) + length(neg)) / ncol(X)
    message(sprintf("cutoff=%s alpha=%s: %d positive, %d negative (%.1f%% of cells)",
                     cutoff, a, length(pos), length(neg), pct * 100))
    if (pct > 0 && pct < cutoff) {
      predicted[pos] <- 1
      predicted[neg] <- -1
      scores <- Coefs
      names(scores) <- colnames(X)
      converged <- TRUE
      used_cutoff <- cutoff
      break
    }
  }
  if (converged) break
}
if (!converged) {
  message("Scissor did not converge to a non-trivial solution at any (cutoff, alpha) combination tried.")
} else if (used_cutoff != 0.2) {
  message(sprintf("NOTE: converged only at relaxed cutoff=%s (0.2 did not converge for any alpha).", used_cutoff))
}

write.csv(data.frame(cell_id = names(predicted), predicted = predicted, score = scores, used_cutoff = used_cutoff),
          out_path, row.names = FALSE)
message(sprintf("Wrote predictions for %d cells to %s", length(predicted), out_path))
