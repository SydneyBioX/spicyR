#' Check that spicy() is calibrated on your data by shuffling the condition
#'
#' Shuffles the condition across patients (or across images, without `subject`), which removes any real difference
#' between the groups while keeping everything else about the study, and runs [spicy()] on each shuffle. If the test
#' is calibrated, about 5% of pairs have p < 0.05 and about 1% p < 0.01 in each shuffle. Pairs share cells, so the
#' percentage moves by a few points from one shuffle to the next; look at the average over shuffles.
#'
#' @param cells The cells, as for [spicy()]: a `SpatialExperiment`, `SingleCellExperiment` or `data.frame` (not a
#'   store).
#' @param condition The column of the condition to shuffle (a factor, character or `Surv` column).
#' @param subject The column of the patient. Each patient keeps one condition, moved as a whole; without it, the
#'   condition is shuffled across images.
#' @param nShuffles The number of shuffles.
#' @param thresholds The p-value thresholds to report.
#' @param ... Other arguments of [spicy()], such as `r`, `imageID`, `cellType` or `covariates`.
#'
#' @return A `data.frame` with one row per shuffle and threshold: `shuffle`, `threshold` and `percent`, the
#'   percentage of pairs whose p-value is below the threshold. `summary()` gives the average percentage over the
#'   shuffles at each threshold, with the expected one, and `plot()` shows each shuffle against the expected
#'   percentage.
#'
#' @examples
#' data("diabetesData")
#' set.seed(1)
#' shuffles <- calibrationCheck(diabetesData, condition = "stage", subject = "case", r = 50, nShuffles = 2)
#' summary(shuffles)
#' plot(shuffles)
#' @export
#' @concept test
calibrationCheck <- function(cells, condition, subject = NULL, nShuffles = 10, thresholds = c(0.05, 0.01), ...) {
  if (.is_store(cells) || .is_class(cells, "SpatialData"))
    stop("calibrationCheck() needs the cells in memory: a SpatialExperiment, SingleCellExperiment or data.frame.",
         call. = FALSE)
  imageID <- list(...)$imageID %||% "imageID"
  unit <- subject %||% imageID
  cd <- if (.is_class(cells, "SummarizedExperiment")) SummarizedExperiment::colData(cells) else cells
  miss <- setdiff(c(condition, unit), names(cd))
  if (length(miss)) stop("columns not found: ", paste(miss, collapse = ", "), call. = FALSE)
  # one condition per patient (or image), shuffled across those that have one; a missing condition stays missing
  u <- as.character(cd[[unit]])
  first <- !duplicated(u)
  units <- u[first]
  labels <- cd[[condition]][first]
  has <- which(!is.na(labels))
  at <- match(u, units)
  shuffles <- lapply(seq_len(nShuffles), function(i) {
    shuffled <- labels
    shuffled[has] <- labels[has][sample.int(length(has))]
    cells[["shuffledCondition"]] <- shuffled[at]
    res <- spicy(cells, condition = "shuffledCondition", subject = subject, ...)
    p <- if (identical(res$method, "cell")) res$cellResults$p_value else topPairs(res, n = NULL)$p.value
    data.frame(shuffle = i, threshold = thresholds,
               percent = 100 * vapply(thresholds, function(t) mean(p < t, na.rm = TRUE), numeric(1)))
  })
  structure(do.call(rbind, shuffles), class = c("spicyCalibration", "data.frame"))
}

#' @rdname calibrationCheck
#' @param object,x The result of `calibrationCheck()`.
#' @export
summary.spicyCalibration <- function(object, ...) {
  object <- as.data.frame(object)
  out <- stats::aggregate(percent ~ threshold, object, mean)
  out <- out[order(-out$threshold), ]
  data.frame(threshold = out$threshold, expected = 100 * out$threshold, percent = out$percent)
}

#' @rdname calibrationCheck
#' @export
plot.spicyCalibration <- function(x, ...) {
  x <- as.data.frame(x)
  expected <- data.frame(threshold = unique(x$threshold))
  ggplot2::ggplot(x, ggplot2::aes(factor(.data$threshold), .data$percent)) +
    ggplot2::geom_jitter(width = 0.1, height = 0, size = 2, alpha = 0.7) +
    ggplot2::geom_point(ggplot2::aes(y = 100 * .data$threshold), data = expected, shape = 95, size = 14,
                        colour = "#b3261e") +
    ggplot2::scale_y_continuous(limits = c(0, max(15, x$percent)), labels = function(v) paste0(v, "%")) +
    ggplot2::labs(x = "p-value threshold", y = "pairs below the threshold",
                  title = paste(length(unique(x$shuffle)), "shuffles of the condition (red: expected)")) +
    ggplot2::theme_classic()
}
