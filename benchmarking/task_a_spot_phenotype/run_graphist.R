#!/usr/bin/env Rscript
# GRAPHIST Stage 1: bulk expression + phenotype -> spatial-network-regularized
# regression -> Graphist(+/-) spot calls. Ported from
# ../../stage1_bulk_regression_R/BRCA-PACSI/STEP1-BRCA-PACSI.R: same APML1 call with
# penalty="Net", but Omega is a spatial adjacency graph (spot/cell k-NN on spatial
# coordinates) instead of an expression-similarity network -- this spatial regularizer
# is the actual difference between GRAPHIST's Stage 1 and vanilla Scissor.
#
# Omega is NOT computed here: it's read from st_omega.csv, precomputed by
# ../simulate_osmfish_starmap.py using generate_adj_mat() -- the exact same Python
# function used everywhere else in GRAPHIST's Stage 2 pipeline -- rather than a second,
# independently-implemented KNN in R, which the original STEP1 scripts never did either
# (they read a precomputed neighbors.csv rather than computing KNN in R).
#
# alpha and cutoff match the original STEP1-BRCA-PACSI.R (alpha=0.03, cutoff=0.2) as the
# first try, extended to a small search only if that exact original setting doesn't
# converge -- the original's single fixed alpha was manually tuned per-dataset, which
# doesn't generalize to new synthetic scenarios, but we keep it as the first, most
# faithful attempt rather than starting from a search.
# Usage: Rscript run_graphist.R <data_dir> <out_csv> [family]

args <- commandArgs(trailingOnly = TRUE)
data_dir <- args[1]
out_path <- args[2]
family <- if (length(args) >= 3) args[3] else "gaussian"

suppressMessages(library(Scissor))
suppressMessages(library(preprocessCore))

bulk <- as.matrix(read.csv(file.path(data_dir, "bulk_expression.csv"), row.names = 1, check.names = FALSE))
phenotype_df <- read.csv(file.path(data_dir, "bulk_phenotype.csv"), row.names = 1)
phenotype <- as.numeric(phenotype_df[, 1])
names(phenotype) <- rownames(phenotype_df)
phenotype <- phenotype[colnames(bulk)]

st <- as.matrix(read.csv(file.path(data_dir, "st_expression.csv"), row.names = 1, check.names = FALSE))
Omega <- as.matrix(read.csv(file.path(data_dir, "st_omega.csv"), row.names = 1, check.names = FALSE))
Omega <- Omega[colnames(st), colnames(st)]

common_genes <- intersect(rownames(bulk), rownames(st))
bulk <- bulk[common_genes, , drop = FALSE]
st <- st[common_genes, , drop = FALSE]

# Quantile-normalize bulk and ST together (matches the original STEP1 scripts).
combined <- cbind(bulk, st)
combined_norm <- normalize.quantiles(as.matrix(combined))
rownames(combined_norm) <- rownames(combined)
colnames(combined_norm) <- colnames(combined)
bulk_norm <- combined_norm[, 1:ncol(bulk), drop = FALSE]
st_norm <- combined_norm[, (ncol(bulk) + 1):ncol(combined_norm), drop = FALSE]

X <- cor(bulk_norm, st_norm)  # bulk samples x ST cells

quality_check <- quantile(X)
message("Bulk-ST correlation five-number summary:")
print(quality_check)

Y <- as.matrix(phenotype)

# Original STEP1-BRCA-PACSI.R value first, then a small search if it doesn't converge.
alpha_seq <- c(0.03, 0.1, 0.3, 0.5)
# Try cutoff=0.2 (the convention used throughout stage1_bulk_regression_R/) first; only
# relax it if nothing in alpha_seq converges under it. Whichever cutoff actually worked is
# logged and written into the output, so relaxed-cutoff results are never silently mixed in
# with 0.2 results.
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
  message("GRAPHIST Stage 1 did not converge to a non-trivial solution at any (cutoff, alpha) combination tried.")
} else if (used_cutoff != 0.2) {
  message(sprintf("NOTE: converged only at relaxed cutoff=%s (0.2 did not converge for any alpha).", used_cutoff))
}

write.csv(data.frame(cell_id = names(predicted), predicted = predicted, score = scores, used_cutoff = used_cutoff),
          out_path, row.names = FALSE)
message(sprintf("Wrote predictions for %d cells to %s", length(predicted), out_path))
