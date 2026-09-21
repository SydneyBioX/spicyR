# Time the R implementation on a cells CSV and write its outputs.
#
# Usage: Rscript benchmarks/run_r.R <spicyR> <cells.csv> <out_dir> <name> <family> <r_or_k> <diagnostics> [tight] [cr2Method]
#
# Prints one line "elapsed_seconds=<s>". Passing "tight" refits iteratively fitted
# pairs to a tight convergence tolerance for output comparison (do not time those
# runs); pass "-" to skip it. cr2Method defaults to "fast".
# Needs binom (diagnostics) and brglm2 (binomial).

args <- commandArgs(trailingOnly = TRUE)
here <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))))
source(file.path(here, "..", "tests", "r_reference", "run_spicyglm.R"))
load_spicyR(args[1])
if (length(args) >= 8 && args[8] == "tight") tighten_iterative_fits()
cr2_method <- if (length(args) >= 9) args[9] else "fast"

cells <- read.csv(args[2])
family <- args[5]
size <- as.numeric(args[6])
elapsed <- run_spicyglm(cells, args[3], args[4], subject = "subject", estimator = "firth",
                        r = if (family == "poisson") size else NULL,
                        k = if (family == "binomial") as.integer(size) else NULL,
                        family = family, diagnostics = as.logical(args[7]), cr2Method = cr2_method)
cat(sprintf("elapsed_seconds=%.3f\n", elapsed))
