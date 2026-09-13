# =============================================================================
# Step-1 check for the observed/expected proportion ratio from getPairwiseProp().
#
# Question: does calcWeights() produce usable (non-NA) weights for the ratio
# statistic v = pi_hat / p0 on REAL data (diabetesData), or does it hit the same
# scale mismatch that broke the arcsine version?
#
# calcWeights() (R/spicy.R:1049) fits  scam( log10(resSq + 1) ~ s(counts) ),
# predicts z1, and floors with  quantile(z1[z1 > 0.1], 0.01).  If nothing clears
# z1 > 0.1 the quantile is NA -> every weight NA -> failed fit -> all-NA p-values.
#
# Run from the spicyR repo root:  Rscript scratch_reproduce_NA_weights.R
# =============================================================================

suppressPackageStartupMessages(devtools::load_all(".", quiet = TRUE))
data("diabetesData")

fmt    <- spicyR:::.format_data(diabetesData, "imageID", "cellType", c("x", "y"), FALSE)
nCells <- table(spicyR:::getImageID(fmt), spicyR:::getCellType(fmt))
from   <- as.character(unique(spicyR:::getCellType(fmt)))
m1     <- rep(from, times = length(from))
m2     <- rep(from, each = length(from))

mL    <- suppressWarnings(suppressMessages(getPairwise(fmt)))       # L(r) - r
mProp <- getPairwiseProp(fmt, k = 15)                               # v = pi_hat / p0

cat(sprintf("dims: L %s   ratio %s\n",
            paste(dim(mL), collapse = "x"), paste(dim(mProp), collapse = "x")))
cat(sprintf("ratio: %.1f%% NA, non-NA range %.3g .. %.3g, median %.3f\n\n",
            100 * mean(is.na(mProp)),
            min(mProp, na.rm = TRUE), max(mProp, na.rm = TRUE),
            median(mProp, na.rm = TRUE)))

# ---- pooled calcWeights() trace (weightsByPair = FALSE, the default path) ----
resSqOf <- function(mat) apply(mat, 2, function(x) {
  if (sum(!is.na(x)) > 1 && sd(x, na.rm = TRUE) > 0) (x - mean(x, na.rm = TRUE))^2
  else rep(NA_real_, length(x))
})

trace <- function(mat, label) {
  cat("----", label, "----\n")
  rS     <- as.vector(resSqOf(mat))
  count1 <- as.vector(nCells[, m1]); count2 <- as.vector(nCells[, m2])
  keep   <- !is.na(rS)
  rSk <- rS[keep]; c1 <- count1[keep]; c2 <- count2[keep]
  cat(sprintf("  pooled n                     : %d\n", length(rSk)))
  cat(sprintf("  residual^2             range : %.3g .. %.3g\n", min(rSk), max(rSk)))
  cat(sprintf("  log10(residual^2 + 1)   range : %.3g .. %.3g\n",
              min(log10(rSk + 1)), max(log10(rSk + 1))))
  wf <- scam::scam(log10(rSk + 1) ~ s(log10(c1 + 1), bs = "mpd") + s(log10(c2 + 1), bs = "mpd"))
  z1 <- suppressWarnings(stats::predict(
    wf, data.frame(c1 = as.numeric(count1), c2 = as.numeric(count2))))
  cat(sprintf("  scam predictions z1     range : %.3g .. %.3g\n", min(z1), max(z1)))
  cat(sprintf("  sum(z1 > 0.1)                 : %d / %d\n", sum(z1 > 0.1), length(z1)))
  fl <- stats::quantile(z1[z1 > 0.1], 0.01, na.rm = TRUE)
  cat(sprintf("  quantile(z1[z1 > 0.1], 0.01)  : %s\n", format(fl)))
  w <- 1 / pmax(z1, fl); w <- w / mean(w, na.rm = TRUE)
  cat(sprintf("  final weights: %d / %d NA\n\n", sum(is.na(w)), length(w)))
}

cat("================ calcWeights() pooled trace on diabetesData ================\n")
trace(mL,    "L(r) - r  (reference: known to work)")
trace(mProp, "v = pi_hat / p0  (the statistic under test)")

# ---- end-to-end spicy() -------------------------------------------------
n_finite <- function(res) {
  pv <- as.matrix(res$p.value)
  sprintf("%d / %d finite   (cols: %s)", sum(is.finite(pv)), length(pv),
          paste(colnames(res$p.value), collapse = ", "))
}

r_w  <- suppressWarnings(suppressMessages(spicy(
  fmt, condition = "stage", subject = "case", alternateResult = mProp, weights = TRUE)))
r_nw <- suppressWarnings(suppressMessages(spicy(
  fmt, condition = "stage", subject = "case", alternateResult = mProp, weights = FALSE)))

cat("================ end-to-end spicy(alternateResult = ratio) ================\n")
cat("weights = TRUE  :", n_finite(r_w),  "\n")
cat("weights = FALSE :", n_finite(r_nw), "\n")
