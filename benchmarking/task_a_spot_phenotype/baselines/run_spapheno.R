#!/usr/bin/env Rscript
# SpaPheno baseline: bulk cell-type-composition + phenotype -> Elastic Net regression,
# applied to each ST cell's SPACE-smoothed local neighborhood cell-type composition
# (single_cell resolution mode -- no deconvolution needed, since osmFISH/STARmap are
# already single-cell resolution).
#
# Only sources the 3 core functions actually needed (Cell_type_neighborhood,
# BuildPhenoModelAutoAlpha, SpatialKNN) directly from the cloned repo, rather than
# installing the full package -- those only depend on FNN and glmnet (already installed),
# avoiding the package's heavier deps (caret, tidyverse, ggplot2, iml, SpaDo) that aren't
# needed for the core prediction pipeline. Scoring uses our own evaluate.py (same as the
# other two methods) rather than SpaPheno's own Precision_Recall_macroF1, so all three
# methods are scored with identical code.
#
# Note on interpretation: for scenarios where the phenotype-defining group label IS the
# finest cell-type label available (osmfish_medium, starmap_*), the bulk composition
# feature is close to a direct encoding of the phenotype itself (see
# simulate_osmfish_starmap.py's celltype field docs) -- SpaPheno may look artificially
# strong there. Flagged in the benchmarking README, not hidden.
#
# Usage: Rscript run_spapheno.R <data_dir> <out_csv> [family] [k_space]

args <- commandArgs(trailingOnly = TRUE)
data_dir <- args[1]
out_path <- args[2]
family <- if (length(args) >= 3) args[3] else "gaussian"
k_space <- if (length(args) >= 4) as.integer(args[4]) else 50

spapheno_dir <- file.path(data_dir, "..", "..", "repos", "SpaPheno", "R")
suppressMessages(library(glmnet))
suppressMessages(library(FNN))
source(file.path(spapheno_dir, "SpatialKNN.R"))
source(file.path(spapheno_dir, "Cell_type_neighborhood.R"))
source(file.path(spapheno_dir, "BuildPhenoModelAutoAlpha.R"))

bulk_composition <- read.csv(file.path(data_dir, "bulk_composition.csv"), row.names = 1, check.names = FALSE)
phenotype_df <- read.csv(file.path(data_dir, "bulk_phenotype.csv"), row.names = 1)
phenotype <- as.numeric(phenotype_df[, 1])
names(phenotype) <- rownames(phenotype_df)
phenotype <- phenotype[rownames(bulk_composition)]

coords <- read.csv(file.path(data_dir, "st_coords.csv"), row.names = 1)
celltype_df <- read.csv(file.path(data_dir, "st_celltype.csv"), row.names = 1)
celltype <- setNames(as.character(celltype_df[, 1]), rownames(celltype_df))
celltype <- celltype[rownames(coords)]

message(sprintf("Computing SPACE (spatial neighborhood cell-type composition, k=%d)...", k_space))
st_composition <- Cell_type_neighborhood(
  sample_information_coordinate = coords,
  resolution = "single_cell",
  sample_information_cellType = celltype,
  k = k_space
)

# Align columns between bulk and ST composition (fill any cell type missing in one side with 0).
all_types <- union(colnames(bulk_composition), colnames(st_composition))
align_cols <- function(mat, cols) {
  missing <- setdiff(cols, colnames(mat))
  for (m in missing) mat[[m]] <- 0
  as.matrix(mat[, cols, drop = FALSE])
}
bulk_mat <- align_cols(bulk_composition, all_types)
st_mat <- align_cols(as.data.frame(st_composition), all_types)

message("Fitting Elastic Net (BuildPhenoModelAutoAlpha)...")
fit <- BuildPhenoModelAutoAlpha(expr = bulk_mat, pheno = phenotype, family = family)
scores <- as.numeric(predict(fit$model, newx = st_mat, s = fit$lambda))
names(scores) <- rownames(st_mat)

# Same selection-fraction convention as run_graphist.R/run_scissor.R: rank by score,
# take the most extreme `cutoff/2` fraction on each side as positive/negative.
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
