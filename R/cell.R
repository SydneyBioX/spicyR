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

## ---- inputs -----------------------------------------------------------------------------------

## Image-level layout shared by every analysis of a data set: cells sorted by image, integer codes
## for images, types and units (patients), the per-image condition, and the C++ dataset.
.cell_context <- function(cells, condition, subject, image_id, cell_type, spatial_coords, ref = NULL,
                          survival = FALSE) {
  x_col <- spatial_coords[1]; y_col <- spatial_coords[2]
  needed <- c(condition, image_id, cell_type, x_col, y_col, subject)
  missing_cols <- setdiff(needed, names(cells))
  if (length(missing_cols)) stop("columns not found in `cells`: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  image_chr <- as.character(cells[[image_id]])
  image_labels <- sort(unique(image_chr), method = "radix")   # C-locale order, as the Python twin
  image_codes <- match(image_chr, image_labels) - 1L
  ord <- order(image_codes, method = "radix")
  df <- cells[ord, , drop = FALSE]; image_codes <- image_codes[ord]
  n_images <- length(image_labels)
  first <- match(seq_len(n_images) - 1L, image_codes)                       # first row of each image

  # units: patients if `subject` maps several images to one patient, else images
  if (is.null(subject) || length(unique(as.character(cells[[subject]]))) == n_images) {
    unit_labels <- image_labels; image_unit <- seq_len(n_images) - 1L
  } else {
    sub <- as.character(df[[subject]])[first]
    unit_labels <- unique(sub); image_unit <- match(sub, unit_labels) - 1L
  }

  ctx <- list(df = df, image_labels = image_labels, image_codes = image_codes, n_images = n_images,
              unit_labels = unit_labels, image_unit = image_unit, first = first)
  if (!is.null(condition) && !survival) {
    lev <- .cell_levels(cells[[condition]])
    if (length(lev) < 2L)
      stop("`condition` needs at least two levels; found ", length(lev), ".", call. = FALSE)
    if (!is.null(ref)) {
      if (!as.character(ref) %in% lev) stop("`ref` is not a level of `condition`.", call. = FALSE)
      lev <- c(as.character(ref), setdiff(lev, as.character(ref)))
    }
    cond <- as.character(df[[condition]])
    if (any(tapply(cond, image_codes, function(z) length(unique(z))) > 1L))
      stop("'", condition, "' must be constant within each image.", call. = FALSE)
    ctx$levels <- lev
    ctx$image_group <- match(cond[first], lev) - 1L
    if (any(tapply(ctx$image_group, image_unit, function(z) length(unique(z))) > 1L))
      stop("each subject must belong to a single '", condition, "' level.", call. = FALSE)
  }
  ctx$type_labels <- unique(as.character(cells[[cell_type]]))
  type_codes <- as.integer(match(as.character(df[[cell_type]]), ctx$type_labels) - 1L)
  ctx$counts <- unclass(table(factor(image_codes, levels = seq_len(n_images) - 1L),
                              factor(type_codes, levels = seq_along(ctx$type_labels) - 1L)))
  offsets <- as.integer(c(0L, cumsum(tabulate(image_codes + 1L, nbins = n_images))))
  ctx$data <- dataset_create(as.numeric(df[[x_col]]), as.numeric(df[[y_col]]), type_codes, offsets,
                             length(ctx$type_labels))
  ctx
}

## Neighbour sums and the label-clustering factor at one radius (or k).
.cell_graph <- function(ctx, pairs, r = NULL, k = NULL, label_clustering = TRUE, window = "convex", n_threads = 1L) {
  knn <- !is.null(k)
  n_types <- length(ctx$type_labels)
  if (knn) {
    dataset_build_knn(ctx$data, as.integer(k), as.integer(n_threads))
    scale <- stats::median(sqrt(k / (pi * rowSums(ctx$counts) / dataset_image_areas(ctx$data, window))))
    if (label_clustering) dataset_build_radius_index(ctx$data, scale)
    h <- 2 * scale
  } else {
    dataset_build_radius_index(ctx$data, r)
    h <- 2 * r
  }
  g <- list(knn = knn, totals = dataset_pair_neighbour_totals(ctx$data, knn, n_types),
            sq = dataset_pair_neighbour_out_sq_totals(ctx$data, knn, n_types))
  codes <- vapply(pairs, function(p) match(p, ctx$type_labels) - 1L, integer(2))
  g$psi <- if (label_clustering) stats_label_clustering(ctx$data, codes[1, ], codes[2, ], ctx$counts, knn, h)
           else matrix(numeric(0), 0, 0)
  g
}

## One pair's image rows, with the unit and group of each row.
.cell_rows <- function(ctx, g, f, t) {
  rows <- stats_excess_image_rows(g$totals, g$sq, ctx$counts, match(f, ctx$type_labels) - 1L,
                                  match(t, ctx$type_labels) - 1L, g$knn, g$psi)
  rows$unit <- ctx$image_unit[rows$img + 1L]
  if (!is.null(ctx$image_group)) rows$group <- ctx$image_group[rows$img + 1L]
  rows
}

## ---- the test: adjusted for abundance and covariates by default ---------------------------

## The log share of the counted type in each image row (the abundance covariate).
.cell_share <- function(ctx, rows, f) {
  i <- rows$img + 1L
  log(pmax(ctx$counts[cbind(i, match(f, ctx$type_labels))], 0.5) / rowSums(ctx$counts)[i])
}

## The design of the main test: one indicator per condition, then the centred log share of the counted
## type (adjust) and the centred covariates. Images with a missing covariate are left out. NULL when there
## is nothing to adjust for (no covariates, and the share is constant or not asked for).
.cell_design <- function(ctx, rows, f, adjust, Z_extra) {
  X <- NULL
  if (adjust) { s <- .cell_share(ctx, rows, f); if (stats::var(s) > 0) X <- cbind(abundance = s) }
  keep <- rep(TRUE, nrow(rows))
  if (!is.null(Z_extra)) { ze <- Z_extra[rows$img + 1L, , drop = FALSE]; keep <- stats::complete.cases(ze); X <- cbind(X, ze) }
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
  r <- stats_excess_test(rows, rows$unit, rows$group, m, frailty, variance)
  out <- list(from = f, to = t, ok = r$ok, reason = r$reason, rows = rows, test = r, unadjusted = r, adjusted_for = "none")
  if (!r$ok) return(out)
  d <- .cell_design(ctx, rows, f, adjust, Z_extra)
  if (is.null(d)) return(out)
  a <- stats_design_tests(d$rows, d$rows$unit, m, d$Z, .cell_contrasts(2L, ncol(d$Z)), .cell_design_tau2(d, 2L, frailty, r$tau2),
                          variance == "hartung_knapp")
  # a design that is not of full rank (e.g. a covariate confounded with the condition): the unadjusted test
  if (!a[[1]]$ok) { out$adjusted_for <- paste0("none (", a[[1]]$reason, ")"); return(out) }
  x <- a[[1]]
  out$test <- list(coef_ref = x$theta[1], coef_comp = x$theta[2], difference = x$estimate, se = x$se, df = x$df,
                   p = x$p, tau2 = x$tau2, influence = x$influence)
  out$effects <- stats::setNames(a[-1], d$extra)
  out$adjusted_for <- d$adjusted_for
  out
}

## More than two conditions: one design with an indicator per level, each level tested against the
## reference (the first level), with the same adjustments as for two conditions.
.cell_pair_test_levels <- function(ctx, rows, f, t, frailty, variance, adjust, Z_extra) {
  G <- length(ctx$levels); m <- length(ctx$unit_labels); hk <- variance == "hartung_knapp"
  per_level <- tapply(rows$unit, factor(rows$group, levels = seq_len(G) - 1L), function(u) length(unique(u)))
  out <- list(from = f, to = t, ok = FALSE, rows = rows, levels = ctx$levels, adjusted_for = "none")
  if (any(is.na(per_level) | per_level < 2L)) { out$reason <- "one_patient_per_group"; return(out) }
  Zg <- outer(rows$group, seq_len(G) - 1L, `==`) * 1
  u <- stats_design_tests(rows, rows$unit, m, Zg, .cell_contrasts(G, G), if (frailty) -1 else 0, hk)
  if (!u[[1]]$ok) { out$reason <- "design_not_full_rank"; return(out) }
  names(u) <- ctx$levels[-1]
  out$ok <- TRUE; out$levels_test <- out$unadjusted_levels <- u
  out$test <- list(coef_ref = u[[1]]$theta[1], tau2 = u[[1]]$tau2)
  d <- .cell_design(ctx, rows, f, adjust, Z_extra)
  if (is.null(d)) return(out)
  a <- stats_design_tests(d$rows, d$rows$unit, m, d$Z, .cell_contrasts(G, ncol(d$Z)),
                          .cell_design_tau2(d, G, frailty, u[[1]]$tau2), hk)
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

.cell_table <- function(fits, ctx, adjusted) {
  if (length(ctx$levels) > 2L) return(.cell_table_levels(fits, ctx, adjusted))
  ok <- vapply(fits, `[[`, logical(1), "ok")
  enames <- .cell_effect_names(fits[ok])
  tab <- do.call(rbind, lapply(fits[ok], function(o) {
    x <- o$test
    row <- data.frame(from = o$from, to = o$to, excess_ref = x$coef_ref, excess_comp = x$coef_comp,
                      excess_difference = x$difference, se = x$se, df = x$df, p_value = x$p, tau2 = x$tau2,
                      stringsAsFactors = FALSE)
    if (!is.null(o$parent)) row$parent <- o$parent
    if (adjusted) {
      row$adjusted_for <- o$adjusted_for
      row <- .cell_add_effects(row, o, enames)
      u <- o$unadjusted
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
  tab <- do.call(rbind, lapply(fits[ok], function(o) {
    do.call(rbind, lapply(names(o$levels_test), function(l) { d <- o$levels_test[[l]]
      row <- data.frame(from = o$from, to = o$to, level = l, excess_ref = d$theta[1],
                        excess_difference = d$estimate, se = d$se, df = d$df, p_value = d$p, tau2 = d$tau2,
                        stringsAsFactors = FALSE)
      if (!is.null(o$parent)) row$parent <- o$parent
      if (adjusted) {
        row$adjusted_for <- o$adjusted_for
        row <- .cell_add_effects(row, o, enames)
        u <- o$unadjusted_levels[[l]]
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
  first <- intersect(c("from", "to", "parent", "level", "r", "excess_ref", "excess_comp", "excess_difference", "se", "df",
                       "p_value", "p_adj", "p_value_best_radius", "tau2", "adjusted_for"), names(tab))
  un <- grep("^unadjusted_", names(tab), value = TRUE)
  tab[, c(first, setdiff(names(tab), c(first, un)), un), drop = FALSE]
}

## Pair labels: from__to, or from__to__parent for Kontextual triples (a table or a fit).
.cell_labels <- function(x) if (!is.null(x$parent)) paste(x$from, x$to, x$parent, sep = "__") else paste(x$from, x$to, sep = "__")

## Per-image excess (O - E) / n of every pair, for plots and bind(): images x pairs.
.cell_image_excess <- function(fits, ctx) {
  out <- lapply(fits, function(o) {
    v <- rep(NA_real_, ctx$n_images)
    if (!is.null(o$rows) && nrow(o$rows)) v[o$rows$img + 1L] <- (o$rows$O - o$rows$E) / o$rows$n
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
