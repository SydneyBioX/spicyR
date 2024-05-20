################################################################################
#
# Generics for spicy
#
################################################################################


#' A table of the significant results from spicy tests
#'
#' @param x The output from spicy.
#' @param coef Which coefficient to list.
#' @param n Extract the top n most significant pairs.
#' @param adj Which p-value adjustment method to use, argument for p.adjust().
#' @param cutoff A p-value threshold to extract significant pairs.
#' @param figures Round to `figures` significant figures.
#'
#' @return A data.frame
#'
#' @examples
#'
#' data(spicyTest)
#' topPairs(spicyTest)
#'
#' @aliases
#' topPairs,SpicyResults-method
#' topPairs
#' @rdname topPairs
#' @export
setGeneric("topPairs", function(x,
                                coef = NULL,
                                n = 10,
                                adj = "fdr",
                                cutoff = NULL,
                                figures = NULL) {
    standardGeneric("topPairs")
})
setMethod("topPairs", "SpicyResults", function(x,
                                               coef = NULL,
                                               n = 10,
                                               adj = "fdr",
                                               cutoff = NULL,
                                               figures = NULL) {
    if (!methods::is(x, "SpicyResults")) stop("x are not results from spicy")

    if (is.null(coef)) coef <- grep("condition", colnames(x$p.value))[1]

    if (methods::is(coef, "character") & !coef %in% colnames(x$p.value)) stop("coef not a column name")
    if (methods::is(coef, "numeric") & !coef %in% seq_len(ncol(x$p.value))) stop("coef not a column name")
    if (length(coef) > 1) warning("coef needs to be length 1, taking first entry.")
    useCondition <- coef[1]
    pval <- x$p.value[[useCondition]]
    adj.pvalue <- stats::p.adjust(pval, adj)

    comp <- x$comparison

    results <-
        data.frame(
            intercept = x$coefficient[, "(Intercept)"],
            coefficient = x$coefficient[, useCondition],
            p.value = pval,
            adj.pvalue = adj.pvalue,
            from = comp$from,
            to = comp$to
        )
    rownames(results) <- rownames(x$coefficient)
    if (length(results$p.value) > 0) {
        results <- results[order(results$p.value), ]

        if (is.null(cutoff) &&
            !is.null(n)) {
            results <- results[seq_len(pmin(n, nrow(results))), ]
        }
        if (is.null(n) &&
            !is.null(cutoff)) {
            results <- results[which(results$adj.pvalue <= cutoff), ]
        }
        if ((!is.null(n)) &&
            (!is.null(cutoff))) {
            results <- results[which(
                results$adj.pvalue <= cutoff & seq_len(nrow(results)) <= n
            ), ]
        }
    }
    if (!is.null(figures) && is.numeric(figures)) {
        results <- results |> dplyr::mutate(
            dplyr::across(is.numeric, signif, digits = figures)
        )
    }
    results
})
