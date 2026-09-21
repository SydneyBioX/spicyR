# Generate reference results from the R implementation (spicyR, gee branch).
#
# Usage: Rscript tests/r_reference/make_fixtures.R /path/to/spicyR
#
# Writes cells.csv and, per case, results_<case>.csv, skipped_<case>.csv and
# (for diagnostic cases) diag_<table>_<case>.csv into this directory.
# Needs spatstat.geom, dplyr, binom and brglm2.

args <- commandArgs(trailingOnly = TRUE)
spicyR <- if (length(args) >= 1) args[1] else "../spicyR"
out_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))))

source(file.path(out_dir, "run_spicyglm.R"))
load_spicyR(spicyR)
tighten_iterative_fits()

# --- synthetic cells: 2 conditions, 14 subjects, 1-3 images each -----------
set.seed(2024)
rows <- list()
img <- 0
for (s in 1:14) {
  cond <- if (s <= 7) "healthy" else "tumour"
  for (j in seq_len(sample(1:3, 1))) {
    img <- img + 1
    n_cells <- sample(150:400, 1)
    x <- runif(n_cells, 0, 500)
    y <- runif(n_cells, 0, 500)
    type <- sample(c("Tcell", "Bcell", "Tumour", "Macro"), n_cells, replace = TRUE,
                   prob = if (cond == "tumour") c(0.2, 0.2, 0.45, 0.15) else c(0.35, 0.3, 0.2, 0.15))
    # attract T cells toward tumour cells in the tumour condition
    if (cond == "tumour") {
      tum <- which(type == "Tumour"); tc <- which(type == "Tcell")
      if (length(tum) > 0 && length(tc) > 0) {
        anchor <- sample(tum, length(tc), replace = TRUE)
        x[tc] <- x[anchor] + rnorm(length(tc), 0, 15)
        y[tc] <- y[anchor] + rnorm(length(tc), 0, 15)
      }
    }
    # a type that only exists in the healthy condition, to exercise skips
    if (cond == "healthy") type[sample(n_cells, 10)] <- "Rare"
    rows[[length(rows) + 1]] <- data.frame(imageID = sprintf("img%02d", img), subject = sprintf("s%02d", s),
                                           condition = cond, cellType = type, x = x, y = y)
  }
}
cells <- do.call(rbind, rows)
write.csv(cells, file.path(out_dir, "cells.csv"), row.names = FALSE)

# --- neighbour ties: integer lattice points with duplicates, spatstat's nnwhich ---
set.seed(7)
tie_x <- sample(0:25, 500, replace = TRUE)
tie_y <- sample(0:25, 500, replace = TRUE)
tie_nn <- matrix(spatstat.geom::nnwhich(ppp(tie_x, tie_y, window = owin(range(tie_x), range(tie_y)), check = FALSE),
                                        k = 1:8), ncol = 8) - 1  # 0-based
write.csv(data.frame(x = tie_x, y = tie_y, nn = tie_nn), file.path(out_dir, "knn_ties.csv"), row.names = FALSE)

run_spicyglm(cells, out_dir, "subject_convex_firth", subject = "subject", window = "convex",
             estimator = "firth", r = 40, diagnostics = TRUE)
run_spicyglm(cells, out_dir, "image_rectangle_firth", subject = NULL, window = "rectangle",
             estimator = "firth", r = 30, diagnostics = TRUE)
run_spicyglm(cells, out_dir, "subject_convex_mle", subject = "subject", window = "convex",
             estimator = "mle", r = 40)
run_spicyglm(cells, out_dir, "binomial_subject_firth", subject = "subject", estimator = "firth",
             family = "binomial", k = 10)
run_spicyglm(cells, out_dir, "binomial_image_mle", subject = NULL, estimator = "mle",
             family = "binomial", k = 6)
run_spicyglm(cells, out_dir, "naive_subject_convex_firth", subject = "subject", window = "convex",
             estimator = "firth", r = 40, cr2Method = "naive")
run_spicyglm(cells, out_dir, "naive_image_rectangle_mle", subject = NULL, window = "rectangle",
             estimator = "mle", r = 30, cr2Method = "naive")
run_spicyglm(cells, out_dir, "naive_binomial_subject_firth", subject = "subject", estimator = "firth",
             family = "binomial", k = 10, cr2Method = "naive")
