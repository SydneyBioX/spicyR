#' The correlation between the tests of a spicy() result
#'
#' The tests of pairs that share a cell type, or that are nested (a cell type and a group containing it), are
#' correlated through the patients they share. `testCorrelation()` estimates the correlation between every two
#' pairs' tests from the one fit, for a correction for multiple testing that uses it in place of independence.
#'
#' Each pair's difference between conditions is a weighted mean over patients: a patient's summary is the weighted
#' mean of its images' effects (`x$pairwiseAssoc`) with the images' weights (`x$imageWeights`), and its weight is
#' the sum of its images' weights times its condition's information (`x$conditionInformation`).
#'
#' `method = "cr2"`: each patient's CR2 influence on the difference, psi = h (summary - its condition's mean) /
#' sqrt(1 - h) with h the patient's share of its condition's weight, negated in the reference condition. Its
#' squares sum to the test's squared standard error, and the cross-products of two pairs' influences,
#' divided by the two standard errors, estimate the correlation of their differences. It is noisy when a few
#' patients carry most of a pair's weight and one condition is small.
#'
#' `method = "permutation"`: the conditions are permuted over the patients `nPermutations` times; for each
#' permutation every test is recomputed in closed form from the patients' summaries and weights (each condition's
#' weighted mean, its CR2 variance and Satterthwaite df, as [spicy()] computes them; the between-patient variance
#' is held at the fit's, and of [spicy()]'s rules for not testing a pair only the need for two patients in each
#' condition is applied), and its p-value is turned into a z with the sign of its (skewness-corrected) t (0 when a condition
#' has fewer than two of the pair's patients). The correlation of
#' the tests' z over the permutations is returned. The permutations come from a fixed seed, so a result always
#' gives the same matrix, in R and in Python.
#'
#' @param x The output of [spicy()] (or [kontextualEngine()]) with two conditions, one radius (or `k`), no
#'   adjustment (no `covariates`, `adjustAbundance = FALSE`) and `variance = "cr2"`.
#' @param method `"cr2"` (the default) or `"permutation"`.
#' @param nPermutations The number of permutations for `method = "permutation"`.
#' @param cores The number of threads for `method = "permutation"`.
#'
#' @return A symmetric matrix with one row and one column per row of `x$cellResults` (the same row names): the
#'   estimated correlation between the pairs' tests. The diagonal is 1, and rows and columns of pairs not tested (or,
#'   by permutation, whose z is the same in every permutation) are `NA`.
#'
#' @examples
#' data("diabetesData")
#' cells <- diabetesData[, diabetesData$stage != "Long-duration"]
#' res <- spicy(cells, condition = "stage", subject = "case", r = 50, from = c("Tc", "Th"), to = c("beta", "alpha"))
#' testCorrelation(res)
#' testCorrelation(res, method = "permutation", nPermutations = 200)
#'
#' @seealso [spicy()]
#' @export
#' @concept test
testCorrelation <- function(x, method = c("cr2", "permutation"), nPermutations = 1000, cores = 1) {
  method <- match.arg(method)
  if (!is(x, "SpicyResults") || !identical(x$method, "cell"))
    stop("`x` must be the output of spicy() with method = \"cell\".", call. = FALSE)
  if (is.null(x$conditionInformation))
    stop("`x` has no conditionInformation: testCorrelation() needs a test of two conditions (not more conditions, a ",
         "continuous condition or survival) from spicyR 1.99.59 or later.", call. = FALSE)
  if (length(x$r) > 1L)
    stop("testCorrelation() needs one radius (or k): with several radii each pair's test combines them.", call. = FALSE)
  tab <- x$cellResults
  if (!is.null(tab$adjusted_for))
    stop("testCorrelation() needs unadjusted tests (no covariates and adjustAbundance = FALSE).", call. = FALSE)
  if (!is.null(x$variance) && x$variance != "cr2")
    stop("testCorrelation() needs variance = \"cr2\".", call. = FALSE)
  if (method == "permutation" && (!is.numeric(nPermutations) || length(nPermutations) != 1L ||
                                  !is.finite(nPermutations) || nPermutations < 2 || nPermutations > 1e7))
    stop("`nPermutations` must be a number from 2 to 1e7.", call. = FALSE)
  if (!is.numeric(cores) || length(cores) != 1L || !is.finite(cores) || cores < 1)
    stop("`cores` must be a positive number of threads.", call. = FALSE)
  labels <- x$comparisons$labels
  unit <- if (is.null(x$subject)) x$imageID else x$subject
  units <- unique(unit)
  tested <- !is.na(tab$p_value)
  y <- do.call(cbind, unname(x$pairwiseAssoc[labels])); w <- do.call(cbind, unname(x$imageWeights[labels]))
  y[, !tested] <- NA_real_; w[, !tested] <- NA_real_
  r <- stats_test_correlation(y, w, match(unit, units) - 1L, as.integer(x$condition) - 1L, length(units),
                              t(x$conditionInformation[labels, , drop = FALSE]),
                              if (method == "cr2") 0L else as.integer(nPermutations), matrix(0L, 0, 0),
                              as.integer(cores))
  # the per-image values must be those of the reported tests
  ok <- abs(r$difference[tested] - tab$excess_difference[tested]) <= 1e-8 * pmax(abs(tab$excess_difference[tested]), tab$se[tested]) &
    abs(r$se[tested] - tab$se[tested]) <= 1e-8 * tab$se[tested]
  if (!all(ok %in% TRUE))
    stop("the per-image effects and weights of `x` do not give its tests (an adjusted test or a variance other ",
         "than CR2?).", call. = FALSE)
  R <- r$R
  R[is.nan(R)] <- NA_real_
  dimnames(R) <- list(rownames(tab), rownames(tab))
  R
}
