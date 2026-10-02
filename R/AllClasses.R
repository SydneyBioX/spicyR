#' The results of spicy()
#'
#' \code{spicy()} returns a \code{SpicyResults} object: a list with the elements below. Use \code{topPairs()},
#' \code{signifPlot()}, \code{spicyBoxPlot()} and \code{bind()} on it, or take elements with \code{$}.
#'
#' Elements for the cell method (\code{method = "cell"}, the default):
#' \describe{
#'   \item{\code{cellResults}}{the full table, one row per pair (and per level with more than two conditions):
#'     \code{excess_ref}, \code{excess_comp}, \code{excess_difference}, \code{se}, \code{df}, \code{p_value},
#'     \code{p_adj}, \code{tau2}; what the test was adjusted for (\code{adjusted_for}) and the effect and
#'     p-value of each adjustment (\code{abundance_effect}, \code{<covariate>_effect}, \code{<covariate>_p_value});
#'     and the unadjusted test (\code{unadjusted_difference}, \code{unadjusted_se}, \code{unadjusted_df},
#'     \code{unadjusted_p_value}, \code{unadjusted_p_adj}). For survival: \code{score_coefficient},
#'     \code{p_value}, \code{p_adj}, \code{hazard_ratio_sd}, \code{hr_p_value} and the unadjusted p-value and
#'     hazard ratio. With several radii: \code{r}, the radius with the strongest evidence, and
#'     \code{p_value_best_radius}.}
#'   \item{\code{radiusResults}}{with several radii, the table at every radius.}
#'   \item{\code{pairwiseAssoc}}{a list with, for every pair, the excess in every image (\code{NA} where the pair
#'     could not be measured).}
#'   \item{\code{imageWeights}}{a list with, for every pair, each image's weight in the test: its share of its
#'     condition's information (sums to 1 within each condition).}
#'   \item{\code{coefficient}, \code{p.value}, \code{se}, \code{statistic}, \code{df}}{matrices with one row per
#'     pair, as for the image method: the reference excess and the differences, and their tests.}
#'   \item{\code{condition}, \code{imageID}, \code{subject}, \code{nCells}, \code{r}, \code{k}}{the condition and
#'     patient of every image, the cell counts, and the settings of the analysis.}
#' }
#' For the image method (\code{method = "image"}) the object holds the matrices \code{coefficient},
#' \code{p.value}, \code{se}, \code{statistic} and \code{df}, and \code{pairwiseAssoc}, as in spicyR 1.x.
#'
#' @seealso \code{\link{spicy}}, \code{\link{topPairs}}, \code{\link{bind}}
#' @name SpicyResults-class
#' @aliases SpicyResults
#' @exportClass SpicyResults
setClass("SpicyResults", contains = "list")
