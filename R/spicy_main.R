#' Test for changes in the co-localisation of cell types between conditions
#'
#' `spicy()` tests, for every ordered pair of cell types `from` → `to`, whether the co-localisation of
#' the two types differs between conditions, or is associated with survival. A pair asks whether `to`
#' cells are placed near `from` cells more than other cells are.
#'
#' **`method = "cell"` (the default).** In each image, the share of `to` cells with at least one `from`
#' cell within radius `r` is compared with its exact expectation q if the `to` cells were a random choice
#' among the cells that are not `from` cells. The effect (`effect = "allocation"`, the default) is the
#' **fraction of `to` cells placed next to (or kept away from) `from` cells**. For a pair that attracts
#' (more `to` cells next to `from` cells than q over all images together) it is (observed share - q) / (1 - q):
#' if a fraction f of the `to` cells were moved next to `from` cells, the effect is f. For a pair that avoids
#' it is (observed share - q) / q: if a fraction f of the `to` cells that would have a `from` cell nearby were
#' moved away, the effect is -f. Either way it does not depend on how many `from` cells there are or how
#' densely they are packed. The side is chosen once per pair from all images, without the conditions, and is
#' reported in the `side` column. `effect = "count"` gives the number of extra `from`
#' cells within `r` of each `to` cell instead; it also reflects how many `from` cells surround a `to` cell
#' (depth of infiltration), but it grows with how densely the `from` cells are packed, so a change in
#' packing alone can appear as a change in co-localisation. Images are combined within patients and
#' patients within conditions by a frailty GEE with a second random effect for the images of a patient (each image
#' is weighted by 1 / (its sampling variance + the between-image variance), so one large image cannot outweigh the
#' patient's other images), and the difference between conditions is tested with a
#' CR2 cluster-robust variance on Satterthwaite degrees of freedom, with **patients (`subject`) as the
#' units**. The difference is adjusted for any `covariates`, and, with `adjustAbundance = TRUE`, for the
#' log share of the `from` type in each image; the unadjusted test is then reported alongside
#' (`unadjusted_*` columns).
#'
#' When nearly every cell has a `from` cell within `r` (q close to 1), there is little room for attraction and
#' an attracting pair's images carry little information; a smaller `r` is more informative.
#'
#' **`method = "image"` (the original spicyR test).** A per-image L-function summary of each pair is
#' compared between conditions with a weighted linear model, or a mixed model when `subject` is given
#' (Canete et al. 2022).
#'
#' @param cells A `data.frame`, `SingleCellExperiment` or `SpatialExperiment` with one row (column) per
#'   cell, or (cell method) the path of a SpatialData or AnnData Zarr store (read with the Rarr package). From a
#'   store, the columns are read a block of rows at a time and the coordinates a batch of whole images at a time,
#'   so the table is never held in memory; the results are those of the same cells given as a `data.frame`.
#'   Column names then refer to the table's `obs`; `imageID` defaults to the table's region key and
#'   `spatialCoords` to its `obsm` array `"spatial"`. A `SpatialData` object (spatialdataR) read from a store
#'   is read from that store in the same way; one built in memory is analysed from its table.
#' @param condition The column of the image-level condition: two groups, or a `survival::Surv` column
#'   for a survival outcome (from a store, a column of `imageInfo`).
#' @param subject The column of the patient (unit) of each image. Images of one patient are combined;
#'   if omitted, every image is its own patient.
#' @param covariates Image- or patient-level columns to adjust for (cell method: added to the design of
#'   the excess, and the effect of each is reported as `<column>_effect` and `<column>_p_value`;
#'   survival: added to the null Cox model).
#' @param imageID,cellType,spatialCoords Column names of the image, cell type and coordinates. For a Zarr
#'   store, `spatialCoords` is either two `obs` columns or the name of one `obsm` array (its first two columns).
#' @param r Radius (or radii) of the neighbourhood, in the units of the coordinates. Cell method: one
#'   radius (default 50), or several to be combined by `combine`. Image method: the radii of the L
#'   function (default 20, 50 and 100).
#' @param from,to Cell types to test (all ordered pairs by default).
#' @param method `"cell"` (spicyR Cell, the default) or `"image"` (the original spicyR test).
#' @param effect Cell method: `"allocation"` (the default), the extra fraction of `to` cells with at least
#'   one `from` cell within `r`; or `"count"`, the number of extra `from` cells within `r` of each `to` cell.
#'   With `k`, "within `r`" means among the cell's `k` nearest neighbours.
#' @param k Cell method: use the `k` nearest neighbours instead of a radius.
#' @param combine Cell method with several radii: `"maxT"` (max-T with the sandwich correlation across
#'   radii; up to 6 radii), `"davies"` (Davies' bound for the same maximum, for many radii) or `"cauchy"`
#'   (Cauchy combination).
#' @param adjustAbundance Cell method: adjust the test for the log share of the `from` type in each image
#'   (default `FALSE`). Its effect is reported as `abundance_effect`. It does not separate more `from`
#'   cells from more densely packed ones, and it removes real effects when the share tracks the condition.
#' @param variance Cell method: `"cr2"` (the default: CR2 on Satterthwaite df) or `"hartung_knapp"` (the
#'   model-based variance, scaled by the Pearson statistic when that exceeds its df and never below CR2, on m - 2
#'   df).
#' @param ref Cell method: the reference level of `condition`.
#' @param cores Number of threads (cell method) or cores (image method).
#' @param tableKey Cell method with a SpatialData store or object: the name of the cell table in `tables/` (by
#'   default the only table, or the one called `"table"`).
#' @param imageInfo Cell method: a table with one row per image, for `condition`, `subject`, `covariates` or a
#'   survival outcome that are not columns of the cell table: a `data.frame` whose `imageID` column, or else its
#'   row names, gives the image IDs, or the name of a table of the SpatialData store or object. Every image of the
#'   cells needs a row, and every row an image. A `survival::Surv` column given as `condition` makes a survival
#'   analysis, also from a store.
#' @param ... Arguments of the image method: `sigma`, `alternateResult`, `minLambda`, `weights`,
#'   `weightsByPair`, `weightFactor`, `weightZThreshold`, `window`, `window.length`, `edgeCorrect`,
#'   `includeZeroCells`, `verbose`, `BPPARAM`. Supplying `alternateResult` selects the image method.
#' @return A `SpicyResults` object. `topPairs()`, `signifPlot()`, `spicyBoxPlot()` and `bind()` work
#'   for both methods. For the cell method, `$cellResults` holds the full table: for the allocation effect,
#'   the pair's `side` (`"attract"` or `"avoid"`, which sets the scale of the effect), the effect in each
#'   condition (`excess_ref`, `excess_comp`, at the average covariates), the difference
#'   (`excess_difference`), its standard error, df, p-value and BH-adjusted p-value, the frailty variance, the
#'   between-image variance (`tau2_image`) and, when the test was adjusted, what for (`adjusted_for`), the effect
#'   and p-value of each adjustment, and the unadjusted test (`unadjusted_difference`, `unadjusted_p_value`, `unadjusted_p_adj`). `$effect`
#'   records which effect was estimated.
#' @references Canete NP et al. (2022). spicyR: spatial analysis of in situ cytometry data in R.
#'   Bioinformatics 38(11), 3099-3105. \doi{10.1093/bioinformatics/btac268}
#'
#'   Bell RM, McCaffrey DF (2002). Bias reduction in standard errors for linear regression with multi-stage
#'   samples. Survey Methodology 28(2), 169-181.
#'
#'   DerSimonian R, Laird N (1986). Meta-analysis in clinical trials. Controlled Clinical Trials 7(3), 177-188.
#'   \doi{10.1016/0197-2456(86)90046-2}
#'
#'   Pustejovsky JE, Tipton E (2018). Small-sample methods for cluster-robust variance estimation and hypothesis
#'   testing in fixed effects models. Journal of Business & Economic Statistics 36(4), 672-683.
#'   \doi{10.1080/07350015.2016.1247004}
#'
#'   Paule RC, Mandel J (1982). Consensus values and weighting factors. Journal of Research of the National
#'   Bureau of Standards 87(5), 377-385. \doi{10.6028/jres.087.022}
#'
#'   Hartung J, Knapp G (2001). A refined method for the meta-analysis of controlled clinical trials with binary
#'   outcome. Statistics in Medicine 20(24), 3875-3889. \doi{10.1002/sim.1009}
#'
#'   Liu Y, Xie J (2020). Cauchy combination test: a powerful test with analytic p-value calculation under
#'   arbitrary dependency structures. Journal of the American Statistical Association 115(529), 393-402.
#'   \doi{10.1080/01621459.2018.1554485}
#' @examples
#' data("diabetesData")
#' # spicyR Cell: patients ("case") are the units
#' # the extra fraction of Th cells, and of beta cells, with a Tc cell within 50 units
#' res <- spicy(diabetesData, condition = "stage", subject = "case", r = 50,
#'              from = "Tc", to = c("Th", "beta"))
#' topPairs(res)
#' res$cellResults
#'
#' # the extra number of Tc cells within 50 units of each Th cell
#' resCount <- spicy(diabetesData, condition = "stage", subject = "case", r = 50,
#'                   from = "Tc", to = "Th", effect = "count")
#'
#' # the original image-level test
#' resImage <- spicy(diabetesData, condition = "stage", subject = "case",
#'                   from = "Tc", to = "Th", method = "image")
#' topPairs(resImage)
#' @aliases spicy
#' @export
spicy <- function(cells,
                  condition,
                  subject = NULL,
                  covariates = NULL,
                  imageID = "imageID",
                  cellType = "cellType",
                  spatialCoords = c("x", "y"),
                  r = NULL,
                  from = NULL,
                  to = NULL,
                  method = c("cell", "image"),
                  effect = c("allocation", "count"),
                  k = NULL,
                  combine = c("maxT", "davies", "cauchy"),
                  adjustAbundance = FALSE,
                  variance = c("cr2", "hartung_knapp"),
                  ref = NULL,
                  cores = 1,
                  tableKey = NULL,
                  imageInfo = NULL,
                  ...) {
  dots <- list(...)
  if (.is_class(cells, "SpatialData")) {   # spatialdataR: its table's store, or the table in memory
    sd <- .spatialdata_input(cells, tableKey, imageInfo)
    if (!is.null(sd$path)) { cells <- sd$path; tableKey <- sd$table } else {
      cells <- sd$cells; imageInfo <- sd$info
      if (missing(imageID) && !is.null(sd$region_key)) imageID <- sd$region_key
    }
  }
  from_store <- .is_store(cells)
  if (missing(method) && (!is.null(dots$alternateResult) || !is.null(dots$sigma))) {
    method <- "image"
    message("`alternateResult` / `sigma` given: using method = \"image\" (the original spicyR test).")
  }
  method <- match.arg(method)
  if (from_store && method == "image")
    stop("a SpatialData store is read by method = \"cell\"; for method = \"image\" give the cells as a data.frame.",
         call. = FALSE)
  if (method == "image" && !missing(effect)) message("`effect` is used by method = \"cell\" only; ignored.")
  if (method == "image") {
    return(do.call(.spicyImage, c(list(cells = cells, condition = condition, subject = subject, covariates = covariates,
                                       imageID = imageID, cellType = cellType, spatialCoords = spatialCoords,
                                       r = r, from = from, to = to, cores = cores), dots)))
  }
  gone <- intersect(names(dots), c("frailty", "labelClustering", "imageHeterogeneity"))
  if (length(gone)) stop("`", gone[1], "` is not an argument of spicy(): the cell method always uses the patient ",
                         "frailty and the image random effect, without label-clustering inflation.", call. = FALSE)
  if (length(dots)) stop("arguments not used by method = \"cell\": ", paste(names(dots), collapse = ", "),
                         ". They belong to method = \"image\".", call. = FALSE)
  combine <- match.arg(combine); variance <- match.arg(variance); effect <- match.arg(effect)
  if (is.null(r) && is.null(k)) r <- 50
  if (!is.null(k) && length(r)) { message("`k` given: using the k nearest neighbours; `r` is ignored."); r <- NULL }

  # image-level columns that are not columns of the cell table come from imageInfo (one row per image)
  need <- unique(c(condition, subject, covariates))
  if (is.character(imageInfo) && length(imageInfo) == 1L) {
    if (!from_store) stop("`imageInfo` names a table of a SpatialData store or object; give a data.frame here.",
                          call. = FALSE)
    imageInfo <- .store_obs_frame(cells, imageInfo)
  }
  if (from_store) {
    # a SpatialData store: the cell table is read a block of rows at a time, the cells a batch of whole images at a time
    st <- .store_table(cells, tableKey)
    if (missing(imageID) && !is.null(st$attrs$region_key)) imageID <- st$attrs$region_key
    if (missing(spatialCoords) && "spatial" %in% st$obsm) spatialCoords <- "spatial"
    info_cols <- intersect(setdiff(need, st$obs_columns), names(imageInfo))
    if (!is.null(subject) && !subject %in% c(st$obs_columns, info_cols)) stop("`subject` column not found.", call. = FALSE)
    miss <- setdiff(covariates, c(st$obs_columns, info_cols))
    if (length(miss)) stop("covariates not found: ", paste(miss, collapse = ", "), call. = FALSE)
    survival <- condition %in% info_cols && inherits(imageInfo[[condition]], "Surv")
    ctx <- .cell_context_store(st, if (!condition %in% info_cols) condition,
                               if (!is.null(subject) && !subject %in% info_cols) subject, setdiff(covariates, info_cols), imageID, cellType, .store_coord_spec(st, spatialCoords),
                               ref = ref)
  } else {
    # the cells stay in input order (no sorted copy); .cell_context() takes the types in the sorted order
    cells <- .format_data(cells, imageID, cellType, spatialCoords, FALSE, ids = FALSE, sort = FALSE)
    info_cols <- intersect(setdiff(need, names(cells)), names(imageInfo))
    if (!is.null(subject) && !subject %in% c(names(cells), info_cols)) stop("`subject` column not found.", call. = FALSE)
    survival <- inherits(if (condition %in% info_cols) imageInfo[[condition]] else cells[[condition]], "Surv")
    if (survival && !condition %in% info_cols) {
      sv <- cells[[condition]]; cells$.time <- sv[, 1]; cells$.event <- sv[, 2]; cells[[condition]] <- NULL
    }
    ctx <- .cell_context(cells, if (!survival && !condition %in% info_cols) condition,
                         if (!is.null(subject) && !subject %in% info_cols) subject, "imageID", "cellType", c("x", "y"),
                         ref = ref, survival = survival, by_image = TRUE)
  }
  if (!is.null(imageInfo))
    ctx <- .cell_image_info(ctx, .image_info_table(imageInfo, ctx$image_labels, imageID), info_cols, condition,
                            subject, ref, survival)
  pheno <- ctx$pheno
  types <- ctx$type_labels
  bad <- setdiff(c(from, to), types)
  if (length(bad)) stop("cell type not found: ", paste(bad, collapse = ", "), call. = FALSE)
  # from -> to: are `to` cells placed next to `from` cells (allocation: the extra fraction of `to` cells with a `from`
  # cell within r; count: extra `from` cells within r of each `to` cell), beyond the `to` cells being a random subset
  # of the cells that are not `from` cells. This is the core's own (counted, centre) order, so pairs pass through.
  pairs <- enumerate_pairs(from, to, types, "binomial")

  Z_extra <- NULL
  if (!is.null(covariates)) {
    miss <- setdiff(covariates, names(pheno))
    if (length(miss)) stop("covariates not found: ", paste(miss, collapse = ", "), call. = FALSE)
    # one row per image, NA where a covariate is missing (those images are left out of the covariate model)
    Z_extra <- stats::model.matrix(stats::reformulate(covariates), stats::model.frame(stats::reformulate(covariates), pheno,
                                   na.action = stats::na.pass))[, -1, drop = FALSE]
    Z_extra <- sweep(Z_extra, 2, colMeans(Z_extra, na.rm = TRUE))
  }

  radii <- if (is.null(k)) sort(unique(r)) else NA
  if (length(radii) > 1L && !survival && length(ctx$levels) > 2L)
    stop("several radii are supported for two conditions; give one `r`.", call. = FALSE)
  if (combine == "maxT" && length(radii) > 6L && !survival)
    stop("combine = \"maxT\" supports up to 6 radii; use combine = \"davies\" for more.", call. = FALSE)
  if (survival) {
    res <- .cell_survival(ctx, pairs, radii, k, pheno, covariates, cores, adjustAbundance, effect)
  } else {
    # several radii: one neighbour search at the largest radius gives every radius's graph
    graphs <- if (length(radii) > 1L) .cell_graphs_radii(ctx, radii, effect) else
      list(.cell_graph(ctx, pairs, r = if (is.na(radii)) NULL else radii, k = k, label_clustering = FALSE,
                       n_threads = cores, effect = effect))
    if (length(ctx$levels) > 2L) {   # one radius: each pair's design in turn
      g <- graphs[[1]]; g$image_heterogeneity <- TRUE
      fits <- lapply(pairs, function(p) .cell_pair_test(ctx, g, p[1], p[2], TRUE, variance, adjustAbundance, Z_extra))
      res <- list(table = .cell_table(fits, ctx, adjustAbundance || !is.null(covariates)), fits = fits)
    } else res <- .cell_pairs(ctx, graphs, pairs, variance, adjustAbundance, Z_extra, radii, combine, cores)
  }
  out <- .cell_results(res, ctx, pheno, condition, subject, survival, radii, k, effect)
  if (!survival) out$variance <- variance
  out
}

## Two conditions: every pair at every radius in one call to the core (stats_pair_tests: per pair and radius the
## image rows, the image random effect, the frailty test and the adjusted design, as .cell_pair_test() does for one,
## and with several radii their combination), then the tables. Several radii: per-radius tables, and one row per pair
## with the combined p-value (the main test and, when adjusted, the unadjusted test); the other columns are those of
## the radius chosen by max-T. The image rows are those of the radius nearest the median.
.cell_pairs <- function(ctx, graphs, pairs, variance, adjust, Z_extra, radii, combine, n_threads = 1L) {
  allocation <- graphs[[1]]$effect == "allocation"
  K <- length(graphs)
  keep <- if (K > 1L) which.min(abs(radii - stats::median(radii))) else 1L
  from <- vapply(pairs, `[`, "", 1L); to <- vapply(pairs, `[`, "", 2L)
  x <- stats_pair_tests(lapply(graphs, `[[`, if (allocation) "any" else "totals"),
                        lapply(graphs, `[[`, if (allocation) "self_expected" else "sq"), allocation, graphs[[1]]$knn,
                        ctx$counts, graphs[[1]]$psi, match(from, ctx$type_labels) - 1L, match(to, ctx$type_labels) - 1L,
                        ctx$image_unit, ctx$image_group, length(ctx$unit_labels), TRUE, TRUE, variance, adjust,
                        if (is.null(Z_extra)) matrix(0, 0, 0) else Z_extra, keep, if (K > 1L) combine else "none",
                        as.integer(n_threads))
  adjusted <- adjust || !is.null(Z_extra)
  tabs <- lapply(seq_len(K), function(k) .cell_pairs_table(x, k, from, to, allocation, adjust, colnames(Z_extra)))
  fits <- lapply(seq_along(pairs), function(i) list(from = from[i], to = to[i], rows = x$rows[[i]],
                                                    unadjusted = list(image_weight = x$image_weight[[i]])))
  if (K == 1L) return(list(table = tabs[[1]], fits = fits))
  long <- do.call(rbind, Map(function(t, rr) if (!is.null(t)) cbind(r = rr, t), tabs, radii))
  # the rows of `long`: radius by radius, the pairs tested there; one row per pair (in order of first appearance),
  # that of its chosen radius
  pk <- lapply(seq_len(K), function(k) which(x$ok[k, ] == 1L))
  at <- matrix(NA_integer_, K, length(pairs)); at[cbind(rep(seq_len(K), lengths(pk)), unlist(pk))] <- seq_along(unlist(pk))
  key <- unique(unlist(pk))
  tab <- long[at[cbind(x$best[key], key)], , drop = FALSE]
  tab$p_value <- x$comb_p[key]; tab$p_value_best_radius <- x$comb_p_best[key]
  if (adjusted) tab$unadjusted_p_value <- x$comb_u_p[key]
  tab$p_adj <- stats::p.adjust(tab$p_value, "BH")
  if (adjusted) tab$unadjusted_p_adj <- stats::p.adjust(tab$unadjusted_p_value, "BH")
  for (col in c("unadjusted_difference", "unadjusted_se", "unadjusted_df")) tab[[col]] <- NULL
  rownames(tab) <- paste(tab$from, tab$to, sep = "__")
  list(table = .cell_order_columns(tab), fits = fits, radius_table = long)
}

## The table of one radius from stats_pair_tests() (the columns and rows .cell_table() gives for the same fits).
.cell_pairs_table <- function(x, k, from, to, allocation, adjust, zn) {
  s <- which(x$ok[k, ] == 1L)
  if (!length(s)) return(NULL)
  at <- function(m) m[k, s]
  tab <- data.frame(from = from[s], to = to[s], excess_ref = at(x$coef_ref), excess_comp = at(x$coef_comp),
                    excess_difference = at(x$difference), se = at(x$se), df = at(x$df), p_value = at(x$p),
                    tau2 = at(x$tau2), stringsAsFactors = FALSE)
  if (allocation) tab$side <- c("avoid", "attract")[(at(x$side) > 0L) + 1L]
  tab$tau2_image <- at(x$tau2_image)
  if (adjust || !is.null(zn)) {
    d <- at(x$design); ab <- at(x$abundance) == 1L
    with_ab <- adjust & (ab | "abundance" %in% zn)
    lab <- if (is.null(zn)) ifelse(with_ab, "abundance", "") else ifelse(with_ab, "abundance+covariates", "covariates")
    tab$adjusted_for <- ifelse(d == 1L, lab, ifelse(d == 2L, paste0("none (", at(x$design_reason), ")"), "none"))
    # the effect of each added column (abundance first, then the covariates), as named in the fits adjusted
    dk <- which(d == 1L)
    enames <- if (length(dk)) unique(c(if (ab[dk[1]]) "abundance", zn, if (any(ab[dk])) "abundance"))
    for (nm in enames) {
      j <- if (nm == "abundance") ifelse(ab, 1L, 1L + match(nm, zn)) else rep(1L + match(nm, zn), length(s))
      e <- d == 1L & !is.na(j)
      i3 <- cbind(j, k, s)[e, , drop = FALSE]
      v <- rep(NA_real_, length(s)); v[e] <- x$effect[i3]; tab[[paste0(nm, "_effect")]] <- v
      v <- rep(NA_real_, length(s)); v[e] <- x$effect_p[i3]; tab[[paste0(nm, "_p_value")]] <- v
    }
    tab$unadjusted_difference <- at(x$u_difference); tab$unadjusted_se <- at(x$u_se); tab$unadjusted_df <- at(x$u_df)
    tab$unadjusted_p_value <- at(x$u_p)
  }
  tab$p_adj <- stats::p.adjust(tab$p_value, "BH")
  if (adjust || !is.null(zn)) tab$unadjusted_p_adj <- stats::p.adjust(tab$unadjusted_p_value, "BH")
  rownames(tab) <- paste(tab$from, tab$to, sep = "__")
  .cell_order_columns(tab)
}

## Survival: score test and the hazard ratio of the shrunken excess, per pair.
.cell_survival <- function(ctx, pairs, radii, k, pheno, covariates, cores, adjust, effect = "allocation") {
  sv <- .cell_survival_setup(ctx, pheno, covariates)
  rr <- radii[1]
  if (length(radii) > 1L) message("survival uses one radius; using r = ", rr, ".")
  g <- .cell_graph(ctx, pairs, r = if (is.na(rr)) NULL else rr, k = k, label_clustering = FALSE, n_threads = cores,
                   effect = effect)
  g$image_heterogeneity <- TRUE
  fits <- lapply(pairs, function(p) .cell_survival_test(ctx, .cell_rows(ctx, g, p[1], p[2]), p[1], p[2], sv, adjust, covariates))
  list(table = .cell_survival_table(fits, adjust, covariates), fits = fits)
}

## The null Cox model (covariates only) and its martingale residuals, one per patient.
.cell_survival_setup <- function(ctx, pheno, covariates) {
  unit_first <- match(seq_along(ctx$unit_labels) - 1L, ctx$image_unit)
  time <- pheno$.time[unit_first]; event <- as.integer(pheno$.event[unit_first])
  if (any(tapply(pheno$.time, ctx$image_unit, function(z) length(unique(z))) > 1L))
    stop("the survival outcome must be constant within each subject.", call. = FALSE)
  W <- if (is.null(covariates)) matrix(0, length(time), 0) else {
    pu <- pheno[unit_first, , drop = FALSE]
    stats::model.matrix(stats::reformulate(covariates), stats::model.frame(stats::reformulate(covariates), pu,
                        na.action = stats::na.pass))[, -1, drop = FALSE] }
  ok_u <- is.finite(time) & !is.na(event) & stats::complete.cases(W)
  null <- stats_cox_fit(time[ok_u], event[ok_u], W[ok_u, , drop = FALSE])
  if (!null$ok) stop("the null Cox model did not converge.", call. = FALSE)
  M <- rep(NA_real_, length(time)); M[ok_u] <- null$martingale
  list(M = M, time = time, event = event)
}

## The survival test of one pair on given image rows.
.cell_survival_test <- function(ctx, rows, f, t, sv, adjust, covariates) {
  m <- length(ctx$unit_labels)
  u <- stats_survival_test(rows, rows$unit, m, sv$M, sv$time, sv$event, numeric(0))
  x <- if (adjust) .cell_share(ctx, rows, f) else numeric(0)
  s <- if (length(x) && stats::var(x) > 0) stats_survival_test(rows, rows$unit, m, sv$M, sv$time, sv$event, x) else u
  list(from = f, to = t, ok = s$ok, reason = s$reason, rows = rows, surv = s, unadjusted_surv = u,
       adjusted_for = paste(c(if (!identical(s, u)) "abundance", if (!is.null(covariates)) "covariates"), collapse = "+"))
}

.cell_survival_table <- function(fits, adjust, covariates) {
  ok <- vapply(fits, `[[`, TRUE, "ok")
  tab <- do.call(rbind, lapply(fits[ok], function(o) { s <- o$surv
    row <- data.frame(from = o$from, to = o$to, score_coefficient = s$score_coef, score_se = s$score_se, score_df = s$score_df,
               p_value = s$score_p, hazard_ratio_sd = s$hr_sd, log_hr_sd = s$log_hr_sd, log_hr_se = s$hr_se,
               hr_p_value = s$hr_p, log_hr_per_unit = s$log_hr_unit, tau2 = s$tau2,
               adjusted_for = if (nzchar(o$adjusted_for)) o$adjusted_for else "none",
               unadjusted_p_value = if (isTRUE(o$unadjusted_surv$ok)) o$unadjusted_surv$score_p else NA_real_,
               unadjusted_hazard_ratio_sd = if (isTRUE(o$unadjusted_surv$ok)) o$unadjusted_surv$hr_sd else NA_real_,
               stringsAsFactors = FALSE)
    if (!is.null(o$parent)) row <- cbind(row[1:2], parent = o$parent, row[-(1:2)])
    .cell_add_side(row, o) }))
  if (!is.null(tab)) {
    tab$p_adj <- stats::p.adjust(tab$p_value, "BH"); tab$unadjusted_p_adj <- stats::p.adjust(tab$unadjusted_p_value, "BH")
    # the unadjusted columns are the test without the abundance adjustment (covariates enter the null Cox model of both)
    if (!adjust) tab[grep("^unadjusted_", names(tab), value = TRUE)] <- NULL
    if (!adjust && is.null(covariates)) tab$adjusted_for <- NULL
    rownames(tab) <- .cell_labels(tab)
  }
  tab
}

## Assemble a SpicyResults object that topPairs(), signifPlot(), spicyBoxPlot() and bind() understand.
.cell_results <- function(res, ctx, pheno, condition, subject, survival, radii, k, effect = "count") {
  tab <- res$table
  num <- vapply(tab, is.double, TRUE); tab[num] <- lapply(tab[num], function(z) { z[is.nan(z)] <- NA_real_; z })
  if (is.null(tab)) stop("no pair could be tested (each condition needs at least two patients with both cell types).",
                         call. = FALSE)
  out <- list(method = "cell", cellResults = tab)
  if (!survival && length(ctx$levels) > 2L) {
    # wide matrices: one column per level against the reference, as the image method's model terms
    key <- .cell_labels(tab); labels <- unique(key)
    wide <- function(col) { d <- data.frame(row.names = labels)
      d[["(Intercept)"]] <- if (col == "excess_difference") tab$excess_ref[match(labels, key)] else NA_real_
      for (l in ctx$levels[-1]) { k <- tab$level == l; d[[paste0("condition", l)]] <- tab[[col]][k][match(labels, key[k])] }
      d }
    out$coefficient <- wide("excess_difference"); out$p.value <- wide("p_value"); out$se <- wide("se"); out$df <- wide("df")
    out$statistic <- out$coefficient; out$statistic[, -1] <- out$coefficient[, -1] / out$se[, -1]; out$statistic[, 1] <- NA
    out$condition <- factor(ctx$levels[ctx$image_group + 1L], levels = ctx$levels)
    tab <- tab[match(labels, key), ]
  } else {
  labels <- rownames(tab)
  term <- if (survival) "condition" else paste0("condition", ctx$levels[2])
  mk <- function(a, b) { d <- data.frame(a, b); names(d) <- c("(Intercept)", term); rownames(d) <- labels; d }
  if (survival) {
    out$coefficient <- mk(NA_real_, tab$log_hr_sd); out$p.value <- mk(NA_real_, tab$p_value)
    out$se <- mk(NA_real_, tab$score_se); out$statistic <- mk(NA_real_, tab$score_coefficient / tab$score_se)
    out$df <- mk(NA_real_, tab$score_df)
    out$survivalOutcome <- survival::Surv(pheno$.time, pheno$.event)
  } else {
    out$coefficient <- mk(tab$excess_ref, tab$excess_difference); out$p.value <- mk(NA_real_, tab$p_value)
    out$se <- mk(NA_real_, tab$se); out$statistic <- mk(NA_real_, tab$excess_difference / tab$se)
    out$df <- mk(NA_real_, tab$df)
    out$condition <- factor(ctx$levels[ctx$image_group + 1L], levels = ctx$levels)
  }
  }
  labels <- .cell_labels(tab)
  if (!is.null(res$radius_table)) out$radiusResults <- res$radius_table
  out$comparisons <- data.frame(from = tab$from, to = tab$to, labels = labels)
  if (!is.null(tab$parent)) { out$comparisons$parent <- tab$parent; out$isKontextual <- TRUE }
  out$pairwiseAssoc <- .cell_image_excess(res$fits, ctx)[labels]
  out$imageWeights <- .cell_image_weight(res$fits, ctx)[labels]
  out$imageIDs <- ctx$image_labels
  out$imageID <- ctx$image_labels
  if (!is.null(subject)) out$subject <- as.character(pheno[[subject]])
  out$nCells <- ctx$n_cells
  out$alternateResult <- FALSE
  out$r <- if (is.null(k)) radii else NULL
  out$k <- k
  out$effect <- effect
  methods::new("SpicyResults", out)
}
