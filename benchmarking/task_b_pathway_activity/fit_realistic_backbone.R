#!/usr/bin/env Rscript
# Fits scDesign3's per-gene negative-binomial marginal model (with a spatial mean
# surface) to REAL lymph node Visium data, then evaluates the fitted (realistic
# baseline mean, realistic dispersion) at a NEW set of synthetic coordinates -- e.g.
# our simulator's 30x30 grid. This gives a genuinely data-grounded count-noise
# backbone (real per-gene dispersion, real spatial smoothness characteristics) that
# our Python simulator can then inject known pathway activity signal on top of, instead
# of the simpler additive-Gaussian-noise model used in sim1-sim5.
#
# Deliberately does NOT use scDesign3's full copula/gene-gene-correlation step
# (fit_copula/simu_new) -- that's the expensive part and models joint gene-gene
# correlation structure, which isn't needed to validate marginal (per-gene, per-spot)
# pathway-activity recovery. fit_marginal + extract_para alone gives realistic MARGINAL
# count characteristics, which is what we actually need here.
#
# Usage: Rscript fit_realistic_backbone.R <lymph_node_dir> <new_coords_csv> <out_dir> [n_genes] [n_cores]

args <- commandArgs(trailingOnly = TRUE)
lymph_dir <- args[1]
new_coords_path <- args[2]
out_dir <- args[3]
n_genes <- if (length(args) >= 4) as.integer(args[4]) else -1  # -1 = all
n_cores <- if (length(args) >= 5) as.integer(args[5]) else 2

suppressMessages({
  library(scDesign3)
  library(SingleCellExperiment)
})

counts <- as.matrix(read.csv(file.path(lymph_dir, "counts.csv"), row.names = 1, check.names = FALSE))
coords <- read.csv(file.path(lymph_dir, "coords.csv"), row.names = 1)
coords <- coords[colnames(counts), ]

if (n_genes > 0 && n_genes < nrow(counts)) {
  set.seed(0)
  keep <- sort(sample(rownames(counts), n_genes))
  counts <- counts[keep, , drop = FALSE]
}
message(sprintf("Fitting %d genes x %d real lymph node spots...", nrow(counts), ncol(counts)))

sce <- SingleCellExperiment(
  assays = list(counts = counts),
  colData = coords
)

t0 <- Sys.time()
data <- construct_data(
  sce, assay_use = "counts", celltype = NULL, pseudotime = NULL,
  spatial = c("spatial1", "spatial2"), other_covariates = NULL, corr_by = "1"
)
marginal_list <- fit_marginal(
  data = data,
  mu_formula = "s(spatial1, spatial2, bs = 'gp', k = 50)",
  sigma_formula = "1",
  family_use = "nb",
  n_cores = n_cores,
  parallelization = "mcmapply",
  trace = FALSE
)
message(sprintf("fit_marginal done in %.1f min", as.numeric(Sys.time() - t0, units = "mins")))

new_coords <- read.csv(new_coords_path, row.names = 1)
colnames(new_coords) <- c("spatial1", "spatial2")

para <- extract_para(
  sce = sce, assay_use = "counts", marginal_list = marginal_list, n_cores = n_cores,
  family_use = "nb", new_covariate = new_coords, data = data$dat
)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
# extract_para()'s mean_mat/sigma_mat are already spots x genes -- no transpose needed
# (an earlier version of this script incorrectly transposed here, producing genes x spots;
# fixed after catching the mismatch via a row/column identity sanity check downstream).
write.csv(para$mean_mat, file.path(out_dir, "realistic_mu.csv"))     # spots x genes
write.csv(para$sigma_mat, file.path(out_dir, "realistic_sigma.csv")) # spots x genes (NB dispersion)
message(sprintf("Wrote realistic (mu, sigma) for %d genes x %d synthetic spots to %s",
                 nrow(counts), nrow(new_coords), out_dir))
