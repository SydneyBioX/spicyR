## spicyR Cell: the R side of the cell-level analysis. Everything numeric happens in the C++ core
## (src/core); these functions only prepare inputs, call it and assemble tables. The Python twin
## (spicyr) mirrors them line for line. Written with AI assistance (Claude, Anthropic), directed by the
## authors; see NEWS.

## Levels present in a condition column: in level order for a factor, else sorted in C-locale order (as the
## Python twin; R's default sort depends on the locale).
.cell_levels <- function(col) {
  if (is.factor(col)) { present <- unique(as.character(col[!is.na(col)])); return(levels(col)[levels(col) %in% present]) }
  sort(unique(as.character(col[!is.na(col)])), method = "radix")
}

## Every ordered (from, to) pair, from-major (both directions are fitted: the excess is directional).
enumerate_pairs <- function(from, to, all_types, family, parent = NULL) {
  if (is.character(from) && length(from) == 1L && is.character(to) && length(to) == 1L)
    return(list(c(from, to)))
  # Kontextual is directional too, and needs `to` in the context
  if (!is.null(parent)) {
    from_types <- if (is.null(from)) all_types else unique(from)
    to_types <- if (is.null(to)) parent else unique(to)
    return(unlist(lapply(from_types, function(f) lapply(to_types, function(t) c(f, t))),
                  recursive = FALSE))
  }
  # the kNN effect is directional (A->B != B->A), so fit every ordered pair
  if (family == "binomial") {
    from_types <- if (is.null(from)) all_types else unique(from)
    to_types <- if (is.null(to)) all_types else unique(to)
    return(unlist(lapply(from_types, function(f) lapply(to_types, function(t) c(f, t))),
                  recursive = FALSE))
  }
  types <- if (!is.null(from) || !is.null(to)) unique(c(from, to)) else all_types
  cross <- if (length(types) >= 2L)
    utils::combn(types, 2L, simplify = FALSE) else list()
  c(cross, lapply(types, function(t) c(t, t)))
}

## Units: patients if `subject` maps several images to one patient, else images. sv: the subject of each image (NULL
## for none); n_sub: the number of subjects over all the cells (as the images' values when constant within images).
.cell_units <- function(ctx, sv, n_sub = length(unique(as.character(sv)))) {
  if (is.null(sv) || n_sub == ctx$n_images) {
    ctx$unit_labels <- ctx$image_labels; ctx$image_unit <- seq_len(ctx$n_images) - 1L
  } else {
    sub <- as.character(sv)
    ctx$unit_labels <- unique(sub); ctx$image_unit <- match(sub, ctx$unit_labels) - 1L
  }
  ctx
}

## The condition of each image (cv), its levels present (lev; `ref` first) and the image groups. check(): the caller's
## check that the condition is constant within images, after the levels are checked.
.cell_condition <- function(ctx, cv, condition, ref = NULL, lev = .cell_levels(cv), check = function() NULL) {
  if (length(lev) < 2L)
    stop("`condition` needs at least two levels; found ", length(lev), ".", call. = FALSE)
  if (!is.null(ref)) {
    if (!as.character(ref) %in% lev) stop("`ref` is not a level of `condition`.", call. = FALSE)
    lev <- c(as.character(ref), setdiff(lev, as.character(ref)))
  }
  check()
  ctx$levels <- lev
  ctx$image_group <- match(as.character(cv), lev) - 1L
  .cell_check_groups(ctx, condition)
}

.cell_check_groups <- function(ctx, condition) {
  if (!is.null(ctx$image_group) && any(tapply(ctx$image_group, ctx$image_unit, function(z) length(unique(z))) > 1L))
    stop("each subject must belong to a single '", condition, "' level.", call. = FALSE)
  ctx
}

## The per-image table (`imageInfo`): one row per image of the analysis, in its order. Image IDs are the `imageID`
## column, or else the row names. Every image needs a row, and every row an image.
.image_info_table <- function(info, image_labels, imageID) {
  if (!is.data.frame(info))
    stop("`imageInfo` must be a data.frame with one row per image (or the name of a table of the store).", call. = FALSE)
  ids <- if (imageID %in% names(info)) as.character(info[[imageID]]) else rownames(info)
  if (anyDuplicated(ids)) stop("`imageInfo` has more than one row for an image.", call. = FALSE)
  unknown <- setdiff(ids, image_labels)
  if (length(unknown)) stop("images in `imageInfo` not found among the cells: ", paste(utils::head(unknown, 5), collapse = ", "),
                            call. = FALSE)
  m <- match(image_labels, ids)
  if (anyNA(m)) stop("images missing from `imageInfo`: ", paste(utils::head(image_labels[is.na(m)], 5), collapse = ", "),
                     call. = FALSE)
  info[m, , drop = FALSE]
}

## Image-level columns from the per-image table `info` (rows in the analysis' order): `cols` are added to ctx$pheno,
## and the condition (a survival outcome: its time and event) and the subject set when they come from it.
.cell_image_info <- function(ctx, info, cols, condition = NULL, subject = NULL, ref = NULL, survival = FALSE) {
  for (nm in cols) ctx$pheno[[nm]] <- info[[nm]]
  if (!is.null(subject) && subject %in% cols) ctx <- .cell_check_groups(.cell_units(ctx, ctx$pheno[[subject]]), condition)
  if (!is.null(condition) && condition %in% cols) {
    if (survival) {
      sv <- ctx$pheno[[condition]]; ctx$pheno[[condition]] <- NULL
      ctx$pheno$.time <- sv[, 1]; ctx$pheno$.event <- sv[, 2]
    } else ctx <- .cell_condition(ctx, ctx$pheno[[condition]], condition, ref)
  }
  ctx
}

## ---- inputs -----------------------------------------------------------------------------------

## Image-level layout shared by every analysis of a data set: integer codes for images (C-locale order of their IDs),
## types and units (patients), the per-image condition, and the cells per image and type. The cells stay in `cells`, in
## input order, with no sorted copy: one pass in the core gives the cells and the first row of every (image, type), and
## the neighbour sums are taken image by image from the columns themselves (.cell_stream()). The whole-data C++ dataset
## (ctx$data, used by Kontextual and the label-clustering factor) is built only with `dataset = TRUE`. Cell types are
## in order of first appearance; by_image = TRUE: with the cells sorted stably by the levels of the image factor, as
## .format_data() sorts them (spicy() passes its cells unsorted), else in row order.
.cell_context <- function(cells, condition, subject, image_id, cell_type, spatial_coords, ref = NULL,
                          survival = FALSE, by_image = FALSE, dataset = FALSE) {
  x_col <- spatial_coords[1]; y_col <- spatial_coords[2]
  needed <- c(condition, image_id, cell_type, x_col, y_col, subject)
  missing_cols <- setdiff(needed, names(cells))
  if (length(missing_cols)) stop("columns not found in `cells`: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  # images and types as factors (their integer codes are all the core reads); a character column is converted once
  as_factor <- function(v) if (is.factor(v)) v else factor(as.character(v))
  fi <- as_factor(cells[[image_id]]); ft <- as_factor(cells[[cell_type]])
  li <- levels(fi); lt <- levels(ft)
  lay <- cell_block_scan(fi, 1L, length(li), ft, 1L, length(lt))   # (image level, type level + NA): cells and first row
  image_labels <- sort(li[rowSums(lay$counts) > 0], method = "radix")   # C-locale order, as the Python twin
  n_images <- length(image_labels)
  image_map <- match(li, image_labels) - 1L; image_map[is.na(image_map)] <- -1L   # image of each image level
  first_level <- apply(lay$first, 1, function(z) if (all(is.na(z))) NA_integer_ else min(z, na.rm = TRUE))
  first <- first_level[match(image_labels, li)]                              # first row of each image in `cells`

  # units: patients if `subject` maps several images to one patient, else images. A subject constant within images
  # takes the same values in the images' first rows as in all the cells.
  ctx <- list(pheno = cells[first, , drop = FALSE], image_labels = image_labels, n_images = n_images, first = first)
  ctx <- .cell_units(ctx, if (!is.null(subject)) cells[[subject]][first],
                     if (!is.null(subject)) { sv <- cells[[subject]]
                       length(unique(as.character(if (cell_constant_within(sv, fi, first_level)) sv[first] else sv))) })
  if (!is.null(condition) && !survival) {
    cv <- cells[[condition]]
    # constant within images (as the same values): the levels present are those of the images' first rows
    constant <- cell_constant_within(cv, fi, first_level)
    ctx <- .cell_condition(ctx, cv[first], condition, ref, lev = .cell_levels(if (constant) cv[first] else cv),
                           check = function() if (!constant && any(tapply(as.character(cv), image_map[as.integer(fi)],
                             function(z) length(unique(z))) > 1L))
                             stop("'", condition, "' must be constant within each image.", call. = FALSE))
  }
  # types present (the last column of the layout is a missing type), by first appearance
  tcol <- which(colSums(lay$counts) > 0)
  key <- if (by_image) apply(lay$first[, tcol, drop = FALSE], 2, function(z) { i <- which(!is.na(z))[1]; c(i, z[i]) })
         else rbind(apply(lay$first[, tcol, drop = FALSE], 2, min, na.rm = TRUE))
  tcol <- tcol[do.call(order, as.data.frame(t(key)))]
  ctx$type_labels <- c(lt, NA_character_)[tcol]
  type_map <- rep(-1L, length(lt) + 1L); type_map[tcol] <- seq_along(tcol) - 1L   # type of each type level (+ NA)
  # cells per image and type, as table(image code, type code) over the cells gave them
  tb <- table(factor(integer(0), levels = seq_len(n_images) - 1L), factor(integer(0), levels = seq_along(tcol) - 1L))
  tb[] <- lay$counts[match(image_labels, li), tcol, drop = FALSE]
  ctx$counts <- unclass(tb)
  ctx$cells <- list(x = cells[[x_col]], y = cells[[y_col]], image = fi, image_map = image_map, type = ft,
                    type_map = type_map)   # the columns themselves, not copies
  if (dataset) {
    ord <- order(image_map[as.integer(fi)], method = "radix")   # cells sorted (stably) by image
    ctx$ord <- ord   # the dataset's rows, as rows of `cells`
    tc <- as.integer(ft)[ord]; tc[is.na(tc)] <- length(lt) + 1L
    ctx$data <- dataset_create(as.numeric(cells[[x_col]][ord]), as.numeric(cells[[y_col]][ord]), type_map[tc],
                               as.integer(c(0L, cumsum(rowSums(ctx$counts)))), length(ctx$type_labels))
  }
  # cells per image and type as table(imageID, cellType) over the cells would give (factor levels, all levels), from
  # the counts already made, so the full columns need not be tabulated again
  lv <- function(x, f) if (is.factor(x)) levels(x) else sort(levels(f))   # levels(f): the unique values
  li <- lv(cells[[image_id]], fi); lt <- lv(cells[[cell_type]], ft)
  nc <- ctx$counts[match(li, image_labels), match(lt, ctx$type_labels), drop = FALSE]; nc[is.na(nc)] <- 0L
  ctx$n_cells <- structure(array(as.integer(nc), dim = dim(nc), dimnames = stats::setNames(list(li, lt), c("", ""))), class = "table")
  ctx
}

## The neighbour totals of every radius (or k), image by image from the cells' own columns (cell_stream_totals): one
## element per radius, list(any, self_expected, self_var) or list(totals, sq, self_var), as the whole-data dataset
## gives them (self_var: the exact self-pair variances, types x images).
## How the allocation self-pair variance's shared-neighbour term is computed (spicy()'s selfVariance and
## selfVarianceCells, kept in ctx$self_var): the core's method code and number of sampled cells. Radius graph only;
## kNN and the count effect are always exact.
.self_var_opts <- function(ctx) {
  s <- if (is.null(ctx$self_var)) list(method = "sampled", cells = 50L) else ctx$self_var
  list(method = match(s$method, c("exact", "geometric", "sampled", "geometricSampled")) - 1L, cells = as.integer(s$cells))
}

.cell_stream <- function(ctx, radii = NULL, k = NULL, allocation = TRUE, n_threads = 1L) {
  if (!is.null(ctx$store)) return(.store_stream(ctx, radii = radii, k = k, allocation = allocation, n_threads = n_threads))
  z <- ctx$cells
  cell_stream_totals(z$x, z$y, z$image, z$image_map, ctx$n_images, z$type, z$type_map, length(ctx$type_labels),
                     as.numeric(if (is.null(radii)) numeric(0) else radii), if (is.null(k)) 0L else as.integer(k),
                     allocation, as.integer(n_threads), .self_var_opts(ctx)$method, .self_var_opts(ctx)$cells)
}

## Neighbour sums and the label-clustering factor at one radius (or k). effect "allocation": for every cell, whether
## it has any cell of each type among its neighbours; "count": how many. Without the label-clustering factor (as
## spicy()) the sums are taken image by image (.cell_stream()); the factor needs the whole-data dataset (ctx$data, from
## .cell_context(dataset = TRUE)).
.cell_graph <- function(ctx, pairs, r = NULL, k = NULL, label_clustering = TRUE, window = "convex", n_threads = 1L,
                        effect = "allocation") {
  knn <- !is.null(k)
  allocation <- effect == "allocation"
  if (!label_clustering) {
    x <- .cell_stream(ctx, radii = r, k = k, allocation = allocation, n_threads = n_threads)[[1]]
    g <- if (allocation) list(knn = knn, effect = effect, any = x$any, self_expected = x$self_expected)
         else list(knn = knn, effect = effect, totals = x$totals, sq = x$sq)
    g$self_var <- x$self_var
    g$psi <- matrix(numeric(0), 0, 0)
    return(g)
  }
  if (is.null(ctx$data)) stop("the label-clustering factor needs .cell_context(dataset = TRUE).", call. = FALSE)
  n_types <- length(ctx$type_labels)
  if (knn) {
    dataset_build_knn(ctx$data, as.integer(k), as.integer(n_threads))
    scale <- stats::median(sqrt(k / (pi * rowSums(ctx$counts) / dataset_image_areas(ctx$data, window))))
    dataset_build_radius_index(ctx$data, scale)
    h <- 2 * scale
  } else {
    dataset_build_radius_index(ctx$data, r)
    h <- 2 * r
  }
  g <- if (allocation) list(knn = knn, effect = effect, any = dataset_pair_neighbour_any_totals(ctx$data, knn, n_types),
                            self_expected = dataset_self_any_expected(ctx$data, knn, n_types))
       else list(knn = knn, effect = effect, totals = dataset_pair_neighbour_totals(ctx$data, knn, n_types),
                 sq = dataset_pair_neighbour_out_sq_totals(ctx$data, knn, n_types))
  g$self_var <- dataset_self_variances(ctx$data, knn, allocation, n_types, .self_var_opts(ctx)$method, .self_var_opts(ctx)$cells)
  # psi of a `to` type is the median over every counted type (not only the requested pairs), so a pair's result does
  # not depend on which other pairs were asked for
  tos <- unique(vapply(pairs, `[`, "", 2L))
  all_pairs <- expand.grid(from = ctx$type_labels, to = tos, stringsAsFactors = FALSE)
  codes <- rbind(match(all_pairs$from, ctx$type_labels), match(all_pairs$to, ctx$type_labels)) - 1L
  g$psi <- stats_label_clustering(ctx$data, codes[1, ], codes[2, ], ctx$counts, knn, h, allocation)
  g
}

## Several radii from one neighbour search (radius graph, no label-clustering inflation), image by image: one graph
## per radius, each identical to .cell_graph(ctx, pairs, r = radii[i], label_clustering = FALSE, effect = effect).
.cell_graphs_radii <- function(ctx, radii, effect = "allocation") {
  allocation <- effect == "allocation"
  m <- .cell_stream(ctx, radii = radii, allocation = allocation)
  lapply(m, function(x) {
    g <- if (allocation) list(knn = FALSE, effect = effect, any = x$any, self_expected = x$self_expected)
         else list(knn = FALSE, effect = effect, totals = x$totals, sq = x$sq)
    g$self_var <- x$self_var
    g$psi <- matrix(numeric(0), 0, 0)
    g
  })
}

## One pair's image rows, with the unit and group of each row.
.cell_rows <- function(ctx, g, f, t) {
  fc <- match(f, ctx$type_labels) - 1L; tc <- match(t, ctx$type_labels) - 1L
  # images with a missing condition give no row: they play no part in the side, tau2_image or the test
  rows <- if (g$effect == "allocation") stats_allocation_image_rows(g$any, g$self_expected, ctx$counts, fc, tc, g$psi, g$self_var,
                                                                    ctx$image_group)
          else stats_excess_image_rows(g$totals, g$sq, ctx$counts, fc, tc, g$knn, g$psi, g$self_var, ctx$image_group)
  rows$unit <- ctx$image_unit[rows$img + 1L]
  if (!is.null(ctx$image_group)) rows$group <- ctx$image_group[rows$img + 1L]
  # allocation: the pair's side, chosen in the core from the images (attraction scales by 1 - q, avoidance by q)
  side <- if (g$effect == "allocation") { if (sum(rows$O - rows$E) < 0) "avoid" else "attract" }
  # images within a patient: each row becomes its effect (O - E) / n with variance v / n^2 + tau2_image
  if (isTRUE(g$image_heterogeneity) && nrow(rows)) {
    h <- stats_image_heterogeneity_rows(rows, rows$unit)
    h$unit <- rows$unit; h$group <- rows$group
    rows <- h
  }
  attr(rows, "side") <- side
  rows
}

## spicy() with soft cell type membership (cellWeights): for each `from` label, one neighbour pass gives the
## weighted random-labelling sums of every `to` label (Dataset::weighted_excess_sums); each pair's image rows then go
## through the side, the image random effect and the test as .cell_pair_test() takes them for cell types. With 0/1
## weights the results are those of the cell types. Self-pairs are not tested.
.cell_weighted <- function(ctx, memb, pairs, r, effect, variance, adjust, Z_extra, survival = FALSE, sv = NULL,
                           covariates = NULL) {
  allocation <- effect == "allocation"
  dataset_build_radius_index(ctx$data, r)
  csr <- .membership_rows(memb, ctx$ord)
  K <- length(memb$labels)
  pairs <- Filter(function(p) p[1] != p[2], pairs)
  if (!length(pairs)) stop("with `cellWeights`, `from` and `to` must give at least one pair of different labels.",
                           call. = FALSE)
  from <- vapply(pairs, `[`, "", 1L); to <- vapply(pairs, `[`, "", 2L)
  fits <- vector("list", length(pairs))
  wctx <- .weighted_context(ctx, memb)
  for (f in unique(from)) {
    jf <- match(f, memb$labels)
    a <- numeric(memb$nCells); sel <- memb$j == jf; a[memb$i[sel]] <- memb$x[sel]
    sums <- dataset_weighted_excess_sums(ctx$data, a[ctx$ord], csr$ptr, csr$lab, csr$w, K, allocation)
    for (k in which(from == f)) {
      t <- to[k]
      rows <- stats_weighted_image_rows(sums, K, match(t, memb$labels) - 1L, allocation, ctx$image_group)
      rows$unit <- ctx$image_unit[rows$img + 1L]
      if (!is.null(ctx$image_group)) rows$group <- ctx$image_group[rows$img + 1L]
      side <- if (allocation) { if (sum(rows$O - rows$E) < 0) "avoid" else "attract" }
      if (nrow(rows)) {
        h <- stats_image_heterogeneity_rows(rows, rows$unit)
        h$unit <- rows$unit; h$group <- rows$group
        rows <- h
      }
      attr(rows, "side") <- side
      fits[[k]] <- if (survival) .cell_survival_test(wctx, rows, f, t, sv, adjust, covariates)
                   else .cell_rows_test(wctx, rows, f, t, TRUE, variance, adjust, Z_extra)
    }
  }
  tab <- if (survival) .cell_survival_table(fits, adjust, covariates)
         else .cell_table(fits, ctx, adjust || !is.null(Z_extra))
  list(table = tab, fits = fits)
}

## Add the allocation side and the between-image variance of a fit to its table row(s).
.cell_add_side <- function(row, o) {
  s <- attr(o$rows, "side")
  if (!is.null(s)) row$side <- s
  h <- attr(o$rows, "tau2_image")
  if (!is.null(h)) row$tau2_image <- h
  row
}

## ---- the test: adjusted for abundance and covariates by default ---------------------------

## The log share of the counted type in each image row (the abundance covariate).
.cell_share <- function(ctx, rows, f) {
  i <- rows$img + 1L
  total <- if (is.null(ctx$cells_per_image)) rowSums(ctx$counts) else ctx$cells_per_image
  log(pmax(ctx$counts[cbind(i, match(f, ctx$type_labels))], 0.5) / total[i])
}

## The context of the tests with soft membership: counts are each label's summed weight per image (for the abundance
## share), over the image's number of cells (labels may overlap).
.weighted_context <- function(ctx, memb) {
  img <- ctx$cells$image_map[as.integer(ctx$cells$image)][memb$i] + 1L
  counts <- matrix(0, ctx$n_images, length(memb$labels))
  s <- rowsum(memb$x, (img - 1L) * length(memb$labels) + memb$j)
  key <- as.integer(rownames(s)) - 1L
  counts[cbind(key %/% length(memb$labels) + 1L, key %% length(memb$labels) + 1L)] <- s[, 1]
  ctx$cells_per_image <- rowSums(ctx$counts)
  ctx$counts <- counts
  ctx$type_labels <- memb$labels
  ctx
}

## The design of the main test: one indicator per condition, then the centred log share of the counted
## type (adjust) and the centred covariates. Images with a missing condition or covariate are left out, as the
## unadjusted test leaves out images with a missing condition. NULL when there is nothing to adjust for (no
## covariates, and the share is constant or not asked for).
.cell_design <- function(ctx, rows, f, adjust, Z_extra) {
  keep <- if (is.null(rows$group)) rep(TRUE, nrow(rows)) else !is.na(rows$group)
  X <- NULL
  if (adjust) { s <- .cell_share(ctx, rows, f); if (stats::var(s[keep]) > 0) X <- cbind(abundance = s) }
  if (!is.null(Z_extra)) { ze <- Z_extra[rows$img + 1L, , drop = FALSE]; keep <- keep & stats::complete.cases(ze); X <- cbind(X, ze) }
  if (is.null(X)) return(NULL)
  rows <- rows[keep, , drop = FALSE]; X <- sweep(X[keep, , drop = FALSE], 2, colMeans(X[keep, , drop = FALSE]))
  G <- length(ctx$levels)
  list(rows = rows, Z = cbind(outer(rows$group, seq_len(G) - 1L, `==`) * 1, X), extra = colnames(X),
       adjusted_for = paste(c(if (adjust && "abundance" %in% colnames(X)) "abundance", if (!is.null(Z_extra)) "covariates"),
                            collapse = "+"))
}

## tau2 of the adjusted design: re-estimated (exact Paule-Mandel) when every added column is constant within
## patients, else held at the unadjusted value (Supplementary Methods, Remark 3).
.cell_design_tau2 <- function(d, G, frailty, tau2) {
  if (!frailty) return(0)
  X <- d$Z[, -seq_len(G), drop = FALSE]
  patient_level <- all(apply(X, 2, function(z) all(tapply(z, d$rows$unit, function(w) length(unique(w)) == 1L))))
  if (patient_level) -1 else tau2
}

## Contrasts of the adjusted design: the level differences, then one row per added column (its effect).
.cell_contrasts <- function(G, q) {
  C <- matrix(0, G - 1L + q - G, q)
  for (l in 2:G) { C[l - 1L, 1] <- -1; C[l - 1L, l] <- 1 }
  if (q > G) for (j in seq_len(q - G)) C[G - 1L + j, G + j] <- 1
  C
}

.cell_pair_test <- function(ctx, g, f, t, frailty, variance, adjust, Z_extra = NULL)
  .cell_rows_test(ctx, .cell_rows(ctx, g, f, t), f, t, frailty, variance, adjust, Z_extra)

## The test of one pair on given image rows (spicy's excess, or the Kontextual excess).
.cell_rows_test <- function(ctx, rows, f, t, frailty, variance, adjust, Z_extra = NULL) {
  if (length(ctx$levels) > 2L) return(.cell_pair_test_levels(ctx, rows, f, t, frailty, variance, adjust, Z_extra))
  m <- length(ctx$unit_labels)
  rules <- .cell_rules(rows)
  r <- stats_excess_test(rows, rows$unit, rows$group, m, frailty, variance, rules$pooled_dispersion)
  out <- list(from = f, to = t, ok = r$ok, reason = r$reason, rows = rows, test = r, unadjusted = r, adjusted_for = "none")
  if (!r$ok) return(out)
  # not tested (as the core's pair_tests): allocation, a condition with no informative patient (every image at -1 or
  # +1, the bounds of the fraction; its CR2 variance is exactly 0); then a condition whose effective patients fall
  # below min(3, m_g / 2), m_g its patients with rows (its mean rests on a few of them)
  skip <- if (r$se > 0 && .cell_no_information(rows, 2L)) "no_information"
          else if (rules$min_effective_units && r$se > 0 && .cell_few_effective(r$n_eff, r$n_units, rules)) "few_effective_units"
  if (!is.null(skip)) { out$test$p <- NaN; out$unadjusted$p <- NaN }
  out$not_tested <- if (!is.null(skip)) skip else if (!(r$se > 0)) "no_variation"
  d <- .cell_design(ctx, rows, f, adjust, Z_extra)
  if (is.null(d)) return(out)
  if (!is.null(skip)) { out$adjusted_for <- paste0("none (", skip, ")"); return(out) }
  # no variation (the unadjusted test is 0 / 0): an adjusted design would only fit rounding noise; not tested
  if (!(r$se > 0)) {
    out$test$p <- NaN; out$adjusted_for <- "none (no_variation)"
    return(out)
  }
  a <- stats_design_tests(d$rows, d$rows$unit, m, d$Z, .cell_contrasts(2L, ncol(d$Z)), .cell_design_tau2(d, 2L, frailty, r$tau2),
                          variance == "hartung_knapp", if (rules$pooled_dispersion) 2L else 0L)
  # a design that is not of full rank (e.g. a covariate confounded with the condition): the unadjusted test
  if (!a[[1]]$ok) { out$adjusted_for <- paste0("none (", a[[1]]$reason, ")"); return(out) }
  x <- a[[1]]
  out$test <- list(coef_ref = x$theta[1], coef_comp = x$theta[2], difference = x$estimate, se = x$se, df = x$df,
                   p = x$p, tau2 = x$tau2, influence = x$influence)
  out$effects <- stats::setNames(a[-1], d$extra)
  out$adjusted_for <- d$adjusted_for
  out
}

## The rules of the condition test for these rows' effect (the core's test_rules: allocation, image rows with a side, or
## the count effect and Kontextual): the pooled dispersion, and the minimum effective number of patients per condition.
.cell_rules <- function(rows) stats_test_rules(!is.null(attr(rows, "side")))

## TRUE when a condition's effective patients n_eff fall below min(3, m_g / 2), m_g its patients with rows (the core's
## few_effective_units).
.cell_few_effective <- function(n_eff, m, rules) any(n_eff < pmin(rules$min_n_eff, as.numeric(m) / 2))

## Allocation (image rows with a side): TRUE when a condition (of G) has no informative patient, one with at least one
## image whose effect (O - E) / n is not exactly -1 or +1, the bounds of the fraction.
.cell_no_information <- function(rows, G) {
  if (is.null(attr(rows, "side")) || is.null(rows$group)) return(FALSE)
  any(stats_informative_units(rows, rows$unit, rows$group, as.integer(G)) == 0L)
}

## More than two conditions: one design with an indicator per level, each level tested against the
## reference (the first level), with the same adjustments as for two conditions.
.cell_pair_test_levels <- function(ctx, rows, f, t, frailty, variance, adjust, Z_extra) {
  G <- length(ctx$levels); m <- length(ctx$unit_labels); hk <- variance == "hartung_knapp"
  rules <- .cell_rules(rows); ng <- if (rules$pooled_dispersion) G else 0L
  per_level <- tapply(rows$unit, factor(rows$group, levels = seq_len(G) - 1L), function(u) length(unique(u)))
  out <- list(from = f, to = t, ok = FALSE, rows = rows, levels = ctx$levels, adjusted_for = "none")
  if (any(is.na(per_level) | per_level < 2L)) { out$reason <- "one_patient_per_group"; return(out) }
  rg <- rows[!is.na(rows$group), , drop = FALSE]   # images with a missing condition are left out
  Zg <- outer(rg$group, seq_len(G) - 1L, `==`) * 1
  u <- stats_design_tests(rg, rg$unit, m, Zg, .cell_contrasts(G, G), if (frailty) -1 else 0, hk, ng)
  # not of full rank, or no variation (the condition means fit the response exactly; the design test checks it, so no
  # adjusted design is fitted then either)
  if (!u[[1]]$ok) { out$reason <- if (u[[1]]$reason == "no_variation") "no_variation" else "design_not_full_rank"; return(out) }
  if (.cell_no_information(rows, G)) { out$reason <- "no_information"; return(out) }
  if (rules$min_effective_units &&
      .cell_few_effective(stats_effective_units(rg, rg$unit, rg$group, G, u[[1]]$tau2), per_level, rules)) {
    out$reason <- "few_effective_units"; return(out)
  }
  names(u) <- ctx$levels[-1]
  out$ok <- TRUE; out$levels_test <- out$unadjusted_levels <- u
  out$test <- list(coef_ref = u[[1]]$theta[1], tau2 = u[[1]]$tau2)
  d <- .cell_design(ctx, rows, f, adjust, Z_extra)
  if (is.null(d)) return(out)
  a <- stats_design_tests(d$rows, d$rows$unit, m, d$Z, .cell_contrasts(G, ncol(d$Z)),
                          .cell_design_tau2(d, G, frailty, u[[1]]$tau2), hk, ng)
  if (!a[[1]]$ok) { out$adjusted_for <- paste0("none (", a[[1]]$reason, ")"); return(out) }
  out$levels_test <- stats::setNames(a[seq_len(G - 1L)], ctx$levels[-1])
  out$effects <- stats::setNames(a[-seq_len(G - 1L)], d$extra)
  out$test <- list(coef_ref = a[[1]]$theta[1], tau2 = a[[1]]$tau2)
  out$adjusted_for <- d$adjusted_for
  out
}

## ---- assembling the results -----------------------------------------------------------------

## The effect of each added column (abundance first, then the covariates), and the unadjusted test.
.cell_effect_names <- function(fits) unique(unlist(lapply(fits, function(o) names(o$effects))))

.cell_add_effects <- function(row, o, enames) {
  for (nm in enames) {
    e <- o$effects[[nm]]; ok <- !is.null(e) && isTRUE(e$ok)
    row[[paste0(nm, "_effect")]] <- if (ok) e$estimate else NA_real_
    row[[paste0(nm, "_p_value")]] <- if (ok) e$p else NA_real_
  }
  row
}

## rbind() of one-row data.frames whose columns may differ (a pair with no image rows has no tau2_image or side):
## the union of the columns, NA where a row lacks one, in order of first appearance.
.cell_rbind_fill <- function(rows) {
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) return(NULL)
  cols <- unique(unlist(lapply(rows, names)))
  do.call(rbind, lapply(rows, function(r) { for (cn in setdiff(cols, names(r))) r[[cn]] <- NA; r[cols] }))
}

.cell_table <- function(fits, ctx, adjusted) {
  if (length(ctx$levels) > 2L) return(.cell_table_levels(fits, ctx, adjusted))
  ok <- vapply(fits, `[[`, logical(1), "ok")
  enames <- .cell_effect_names(fits[ok])
  # every requested pair, tested or not (not_tested: why not; NA when it was)
  tab <- .cell_rbind_fill(lapply(fits, function(o) {
    x <- if (isTRUE(o$ok)) o$test else list(coef_ref = NA_real_, coef_comp = NA_real_, difference = NA_real_,
                                             se = NA_real_, df = NA_real_, p = NA_real_, tau2 = NA_real_)
    row <- data.frame(from = o$from, to = o$to, excess_ref = x$coef_ref, excess_comp = x$coef_comp,
                      excess_difference = x$difference, se = x$se, df = x$df, p_value = x$p, tau2 = x$tau2,
                      stringsAsFactors = FALSE)
    if (!is.null(o$parent)) row$parent <- o$parent
    row <- .cell_add_side(row, o)
    if (!isTRUE(o$ok) && !is.null(row$side)) row$side <- NA_character_
    row$not_tested <- if (!isTRUE(o$ok)) (if (length(o$reason) && nzchar(o$reason)) o$reason else "not_tested")
                      else if (!is.null(o$not_tested)) o$not_tested else NA_character_
    if (adjusted) {
      row$adjusted_for <- o$adjusted_for
      row <- .cell_add_effects(row, o, enames)
      u <- if (isTRUE(o$ok)) o$unadjusted else list(difference = NA_real_, se = NA_real_, df = NA_real_, p = NA_real_)
      row$unadjusted_difference <- u$difference; row$unadjusted_se <- u$se; row$unadjusted_df <- u$df
      row$unadjusted_p_value <- u$p
    }
    row
  }))
  if (is.null(tab)) return(NULL)
  tab$p_adj <- stats::p.adjust(tab$p_value, "BH")
  if (adjusted) tab$unadjusted_p_adj <- stats::p.adjust(tab$unadjusted_p_value, "BH")
  rownames(tab) <- .cell_labels(tab)
  .cell_order_columns(tab)
}

## More than two conditions: one row per pair and level (contrast with the reference level).
.cell_table_levels <- function(fits, ctx, adjusted) {
  ok <- vapply(fits, `[[`, TRUE, "ok")
  enames <- .cell_effect_names(fits[ok])
  # every requested pair and level, tested or not (not_tested: why not; NA when it was)
  na_test <- list(theta = NA_real_, estimate = NA_real_, se = NA_real_, df = NA_real_, p = NA_real_, tau2 = NA_real_)
  tab <- .cell_rbind_fill(lapply(fits, function(o) {
    lv <- if (isTRUE(o$ok)) names(o$levels_test) else ctx$levels[-1]
    .cell_rbind_fill(lapply(lv, function(l) { d <- if (isTRUE(o$ok)) o$levels_test[[l]] else na_test
      row <- data.frame(from = o$from, to = o$to, level = l, excess_ref = d$theta[1],
                        excess_difference = d$estimate, se = d$se, df = d$df, p_value = d$p, tau2 = d$tau2,
                        stringsAsFactors = FALSE)
      if (!is.null(o$parent)) row$parent <- o$parent
      row <- .cell_add_side(row, o)
      if (!isTRUE(o$ok) && !is.null(row$side)) row$side <- NA_character_
      row$not_tested <- if (isTRUE(o$ok)) NA_character_ else if (length(o$reason)) o$reason else "not_tested"
      if (adjusted) {
        row$adjusted_for <- o$adjusted_for
        row <- .cell_add_effects(row, o, enames)
        u <- if (isTRUE(o$ok)) o$unadjusted_levels[[l]] else na_test
        row$unadjusted_difference <- u$estimate; row$unadjusted_se <- u$se; row$unadjusted_df <- u$df
        row$unadjusted_p_value <- u$p
      }
      row })) }))
  if (is.null(tab)) return(NULL)
  for (l in unique(tab$level)) { k <- tab$level == l
    tab$p_adj[k] <- stats::p.adjust(tab$p_value[k], "BH")
    if (adjusted) tab$unadjusted_p_adj[k] <- stats::p.adjust(tab$unadjusted_p_value[k], "BH") }
  rownames(tab) <- paste(.cell_labels(tab), tab$level, sep = "__")
  .cell_order_columns(tab)
}

## The main test first, then what it was adjusted for and the effects, then the unadjusted test.
.cell_order_columns <- function(tab) {
  first <- intersect(c("from", "to", "parent", "level", "r", "side", "excess_ref", "excess_comp", "excess_difference",
                       "se", "df", "p_value", "p_adj", "p_value_best_radius", "tau2", "tau2_image", "not_tested",
                       "adjusted_for"), names(tab))
  un <- grep("^unadjusted_", names(tab), value = TRUE)
  tab[, c(first, setdiff(names(tab), c(first, un)), un), drop = FALSE]
}

## Pair labels: from__to, or from__to__parent for Kontextual triples (a table or a fit).
.cell_labels <- function(x) if (!is.null(x$parent)) paste(x$from, x$to, x$parent, sep = "__") else paste(x$from, x$to, sep = "__")

## Per-image excess (O - E) / n of every pair (allocation: the extra fraction of `to` cells next to `from`; count: the
## extra `from` cells per `to` cell), for plots and bind(): images x pairs.
.cell_image_excess <- function(fits, ctx) {
  out <- lapply(fits, function(o) {
    v <- rep(NA_real_, ctx$n_images)
    if (length(o$rows$O)) v[o$rows$img + 1L] <- (o$rows$O - o$rows$E) / o$rows$n
    v })
  names(out) <- vapply(fits, .cell_labels, "")
  out
}

## Per-image weight of every pair: the image's share of its condition's information in the frailty model
## (sums to 1 within each condition); NA where the pair was not tested.
.cell_image_weight <- function(fits, ctx) {
  out <- lapply(fits, function(o) {
    v <- rep(NA_real_, ctx$n_images)
    w <- o$unadjusted$image_weight
    if (!is.null(w) && length(w)) v[o$rows$img + 1L] <- w
    v })
  names(out) <- vapply(fits, .cell_labels, "")
  out
}
