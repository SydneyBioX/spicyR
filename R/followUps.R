#' The follow-up comparisons of two levels, for more than two conditions
#'
#' With more than two conditions (cell method), [spicy()] tests each pair first with the overall test (are all the
#' condition means equal? `overall_p_value`, BH-adjusted over the pairs in `overall_p_adj`). For the pairs that test
#' calls (`overall_p_adj` at most `alpha`, R of the m pairs with an overall test), every two levels are then compared:
#' each comparison is the two-condition test of [spicy()] on the patients of those two levels only (their own
#' between-image and between-patient variances, the same rules for not testing, and the same adjustments). A
#' comparison is tested whenever its two levels pass the rules, even when another level of the pair is too sparse.
#' Within a pair the comparisons are Bonferroni-corrected at the Benjamini-Bogomolov level alpha R / m: a comparison
#' is `significant` when its p-value is at most `p_threshold` = alpha R / (m K), K the comparisons tested in the pair.
#' This controls the average, over the pairs called, of the chance of any false comparison within a pair (Benjamini
#' and Bogomolov 2014). Pairs the overall test does not call have no follow-ups.
#'
#' @param x The output from [spicy()] with more than two conditions.
#' @param n The number of comparisons to list (`NULL`, the default, lists every one).
#' @param from,to Only the comparisons of pairs with these `from` and `to` cell types.
#' @param significant `TRUE`: only the comparisons significant at the Benjamini-Bogomolov level.
#' @param figures Round to `figures` significant figures.
#'
#' @return A data.frame with one row per called pair and comparison, ordered by the pair's overall p-value and then
#'   the comparison's: `from`, `to`, the two levels (`level_ref`, `level_comp`), the effect in each
#'   (`excess_ref`, `excess_comp`), their difference (`excess_difference`, comp minus ref), its standard error, df and
#'   p-value, `p_threshold`, `significant`, the pair's `overall_p_value` and `overall_p_adj`, the number of comparisons
#'   tested in the pair (`n_comparisons`), the frailty variance of the comparison (`tau2`) and `not_tested` (why a
#'   comparison was not tested; `NA` when it was). For the allocation effect, the pair's `side` and the between-image
#'   variance of the comparison (`tau2_image`); with adjustment, `adjusted_for` and the unadjusted difference and
#'   p-value. Row names are `"from__to__comp__ref"`. The same table is `x$followUpResults`.
#'
#' @references Benjamini Y, Bogomolov M (2014). Selective inference on multiple families of hypotheses. Journal of
#'   the Royal Statistical Society B 76(1), 297-318. \doi{10.1111/rssb.12028}
#'
#' @examples
#' data("diabetesData")
#' res <- spicy(diabetesData, condition = "stage", subject = "case", r = 50, from = "Tc", to = c("Th", "beta"))
#' topFollowUps(res)
#'
#' @seealso [spicy()], [topPairs()]
#' @export
#' @concept test
topFollowUps <- function(x, n = NULL, from = NULL, to = NULL, significant = FALSE, figures = NULL) {
    if (!identical(x$method, "cell") || is.null(x$cellResults$overall_p_value))
        stop("follow-up comparisons are made by spicy() (cell method) with more than two conditions.", call. = FALSE)
    fu <- x$followUpResults
    if (is.null(fu)) stop("these results have no follow-up comparisons (made by an older spicyR).", call. = FALSE)
    if (!is.null(from)) fu <- fu[fu$from %in% from, , drop = FALSE]
    if (!is.null(to)) fu <- fu[fu$to %in% to, , drop = FALSE]
    if (significant) fu <- fu[which(fu$significant), , drop = FALSE]
    fu <- fu[order(fu$overall_p_value, fu$p_value), , drop = FALSE]
    if (!is.null(n)) fu <- fu[seq_len(min(n, nrow(fu))), , drop = FALSE]
    if (!is.null(figures) && is.numeric(figures)) {
        isNum <- vapply(fu, is.numeric, logical(1))
        fu[isNum] <- lapply(fu[isNum], signif, digits = figures)
    }
    fu
}
