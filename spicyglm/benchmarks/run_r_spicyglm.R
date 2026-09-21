# Time spicyR's spicyGLM() on a cells CSV.
#
# Usage: Rscript benchmarks/run_r_spicyglm.R <spicyR> <cells.csv> <family> <r_or_k> <diagnostics> <cores> <repeats>
#
# Loads spicyR with devtools::load_all() and converts the CSV to a
# SpatialExperiment before timing. Prints "elapsed_seconds=<t1>,<t2>,..." with
# one fitting time per repeat.

args <- commandArgs(trailingOnly = TRUE)
suppressMessages(devtools::load_all(args[1], quiet = TRUE))
suppressPackageStartupMessages(library(SpatialExperiment))

cells <- read.csv(args[2])
spe <- SpatialExperiment(
  colData = cells[, c("imageID", "subject", "condition", "cellType")],
  spatialCoords = as.matrix(cells[, c("x", "y")])
)
family <- args[3]
size <- as.numeric(args[4])
fit_args <- list(cells = spe, condition = "condition", subject = "subject", imageID = "imageID",
                 cellType = "cellType", family = family, cores = as.integer(args[6]),
                 computeDiagnostics = as.logical(args[5]))
if (family == "poisson") fit_args$r <- size else fit_args$k <- as.integer(size)

times <- numeric(0)
for (i in seq_len(as.integer(args[7]))) {
  start <- Sys.time()
  invisible(capture.output(suppressWarnings(suppressMessages(do.call(spicyGLM, fit_args)))))
  times <- c(times, as.numeric(difftime(Sys.time(), start, units = "secs")))
}
cat(sprintf("elapsed_seconds=%s\n", paste(sprintf("%.4f", times), collapse = ",")))
