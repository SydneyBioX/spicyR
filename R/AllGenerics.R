################################################################################
#
# Generics for spicy
#
################################################################################


#' The most significant pairs, or the results of chosen pairs
#'
#' Lists the pairs tested by [spicy()], most significant first. Give `from` and/or `to` to see chosen pairs, and
#' `n = NULL` to see every pair.
#'
#' For the cell method (the default of [spicy()]) the table has one row per pair (and per level with more than two
#' conditions), with the columns of the test: `from`, `to`, `side`, `excess_ref` and `excess_comp` (the average
#' effect in the reference and the comparison group), `excess_difference`, `p_value` and `p_adj`; with several
#' radii, `r`, the radius with the strongest evidence; with covariates, the effect and p-value of each
#' (`<covariate>_effect`, `<covariate>_p_value`); for survival, `hazard_ratio_sd`. `full = TRUE` adds the
#' standard error, degrees of freedom, the variance between patients (`tau2`), what the test was adjusted for and
#' the unadjusted test. For the image method the columns are `intercept`, `coefficient`, `p.value`,
#' `adj.pvalue`, `from` and `to`, as in spicyR 1.x.
#'
#' @param x The output from [spicy()].
#' @param coef Which coefficient to list: a column of `x$p.value`, or with the cell method and more than two
#'   conditions, a level of the condition.
#' @param n The number of pairs to list. `NULL` lists every pair, as does leaving `n` out when `cutoff` is given.
#' @param adj Which p-value adjustment method to use, argument for p.adjust().
#' @param cutoff List only the pairs whose adjusted p-value is at most `cutoff`.
#' @param figures Round to `figures` significant figures.
#' @param from,to Only the pairs with these `from` and `to` cell types.
#' @param full Cell method: every column of the test, not only the main ones.
#'
#' @return A data.frame, one row per pair, with the pair (`"from__to"`) as row names.
#'
#' @examples
#'
#' data(spicyTest)
#' topPairs(spicyTest)
#'
#' data("diabetesData")
#' res <- spicy(diabetesData, condition = "stage", subject = "case", r = 50, from = "Tc", to = c("Th", "beta"))
#' topPairs(res, from = "Tc", to = "Th")
#'
#' @aliases topPairs,SpicyResults-method topPairs
#' @rdname topPairs
#' @export
#' @concept test
setGeneric("topPairs", function(x,
                                coef = NULL,
                                n = 10,
                                adj = "fdr",
                                cutoff = NULL,
                                figures = NULL,
                                from = NULL,
                                to = NULL,
                                full = FALSE) {
    standardGeneric("topPairs")
})
setMethod("topPairs", "SpicyResults", function(x,
                                               coef = NULL,
                                               n = 10,
                                               adj = "fdr",
                                               cutoff = NULL,
                                               figures = NULL,
                                               from = NULL,
                                               to = NULL,
                                               full = FALSE) {
    if (!is.null(cutoff) && missing(n)) n <- NULL
    if (identical(x$method, "cell")) {
        results <- .topCellPairs(x, coef, adj, full)
        pval <- "p_value"
        padj <- "p_adj"
    } else {
        results <- .topImagePairs(x, coef, adj)
        pval <- "p.value"
        padj <- "adj.pvalue"
    }
    results <- .pickPairs(results, from, to)
    results <- results[order(results[[pval]]), , drop = FALSE]
    if (!is.null(cutoff)) results <- results[which(results[[padj]] <= cutoff), , drop = FALSE]
    if (!is.null(n)) results <- results[seq_len(min(n, nrow(results))), , drop = FALSE]
    if (!is.null(figures) && is.numeric(figures)) {
        isNum <- vapply(results, is.numeric, logical(1))
        results[isNum] <- lapply(results[isNum], signif, digits = figures)
    }
    results
})

# the rows of the chosen from and to cell types
.pickPairs <- function(results, from, to) {
    want <- list(from = from, to = to)
    for (side in names(want)[!vapply(want, is.null, TRUE)]) {
        miss <- setdiff(want[[side]], results[[side]])
        if (length(miss)) stop("no pair with ", side, " = ", paste(dQuote(miss, FALSE), collapse = ", "),
                               " in the results.", call. = FALSE)
        results <- results[results[[side]] %in% want[[side]], , drop = FALSE]
    }
    results
}

# the cell method's table: the main columns, or with full = TRUE all of them
.topCellPairs <- function(x, coef, adj, full) {
    tab <- x$cellResults
    if (!is.null(coef) && !is.null(tab$level)) tab <- tab[tab$level == sub("^condition", "", coef), , drop = FALSE]
    if (!adj %in% c("fdr", "BH")) tab$p_adj <- stats::p.adjust(tab$p_value, adj)
    if (full) return(tab)
    effects <- grep("_effect$", names(tab), value = TRUE)
    covariate <- c(rbind(effects, sub("_effect$", "_p_value", effects)))
    main <- c("from", "to", "parent", "level", "r", "side", "excess_ref", "excess_comp", "excess_difference",
              "hazard_ratio_sd", "p_value", "p_adj", covariate)
    tab[intersect(main, names(tab))]
}

# the image method's table, as in spicyR 1.x
.topImagePairs <- function(x, coef, adj) {
    if (is.null(x$p.value) && !is.null(x$survivalResults)) {   # survival: one coefficient, the log hazard ratio
        s <- x$survivalResults
        pair <- do.call(rbind, strsplit(s$test, "__", fixed = TRUE))
        return(data.frame(coefficient = s$coef, se = s$se.coef, p.value = s$p.value,
                          adj.pvalue = stats::p.adjust(s$p.value, adj), from = pair[, 1], to = pair[, 2],
                          row.names = s$test))
    }
    if (is.null(coef)) coef <- grep("condition", colnames(x$p.value))[1]
    if (methods::is(coef, "character") & !coef %in% colnames(x$p.value)) stop("coef not a column name")
    if (methods::is(coef, "numeric") & !coef %in% seq_len(ncol(x$p.value))) stop("coef not a column name")
    if (length(coef) > 1) warning("coef needs to be length 1, taking first entry.")
    useCondition <- coef[1]
    pval <- x$p.value[[useCondition]]
    comp <- x$comparisons
    results <- data.frame(
        intercept = x$coefficient[, "(Intercept)"],
        coefficient = x$coefficient[, useCondition],
        p.value = pval,
        adj.pvalue = stats::p.adjust(pval, adj),
        from = comp$from,
        to = comp$to
    )
    if (!is.null(comp$parent)) results$parent <- comp$parent
    rownames(results) <- rownames(x$coefficient)
    results
}
