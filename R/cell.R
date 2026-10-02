## spicyR Cell: the R side of the cell-level analysis. Everything numeric happens in the C++ core
## (src/core); these functions only prepare inputs, call it and assemble tables. The Python twin
## (spicyr) mirrors them line for line.

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

## ---- the test, with the availability adjustment and covariates ---------------------------------

.cell_pair_test <- function(ctx, g, f, t, frailty, variance, availability, Z_extra = NULL) {
  rows <- .cell_rows(ctx, g, f, t)
  if (length(ctx$levels) > 2L) return(.cell_pair_test_levels(ctx, rows, f, t, frailty, availability, Z_extra))
  r <- stats_excess_test(rows, rows$unit, rows$group, length(ctx$unit_labels), frailty, variance)
  out <- list(from = f, to = t, ok = r$ok, reason = r$reason, rows = rows, test = r)
  if (!r$ok) return(out)
  if (!is.null(Z_extra)) {
    # covariates: the group difference adjusted for image-level covariates; tau2 re-estimated when
    # every covariate is constant within patients (exact Paule-Mandel), else held (new_methods.pdf, Remark 3).
    # Images with a missing covariate are left out of this model.
    ze <- Z_extra[rows$img + 1L, , drop = FALSE]; cc <- stats::complete.cases(ze)
    rc <- rows[cc, , drop = FALSE]; ze <- ze[cc, , drop = FALSE]
    Z <- cbind(rc$group == 0, rc$group == 1, ze)
    patient_level <- all(apply(ze, 2, function(z) all(tapply(z, rc$unit, function(w) length(unique(w)) == 1L))))
    out$covariate <- stats_design_test(rc, rc$unit, length(ctx$unit_labels), Z, c(-1, 1, rep(0, ncol(Z_extra))),
                                       if (patient_level) -1 else r$tau2)
  }
  if (availability) {
    i <- rows$img + 1L
    share <- log(pmax(ctx$counts[cbind(i, match(f, ctx$type_labels))], 0.5) / rowSums(ctx$counts)[i])
    out$availability <- if (stats::var(share) > 0)
      stats_availability_test(rows, rows$unit, rows$group, length(ctx$unit_labels), share, r$tau2)
    else list(ok = FALSE, reason = "share_constant")
  }
  out
}

## More than two conditions: one design with an indicator per level, each level tested against the
## reference (the first level). The design is patient-level, so Paule-Mandel is exact; tau2 is
## estimated once and shared by the level contrasts (CR2 on Satterthwaite df for each).
.cell_pair_test_levels <- function(ctx, rows, f, t, frailty, availability, Z_extra) {
  G <- length(ctx$levels); m <- length(ctx$unit_labels)
  per_level <- tapply(rows$unit, factor(rows$group, levels = seq_len(G) - 1L), function(u) length(unique(u)))
  out <- list(from = f, to = t, ok = FALSE, rows = rows, levels = ctx$levels)
  if (any(is.na(per_level) | per_level < 2L)) { out$reason <- "one_patient_per_group"; return(out) }
  Zg <- outer(rows$group, seq_len(G) - 1L, `==`) * 1
  fit <- function(Z, tau2) {
    q <- ncol(Z); tau <- tau2; res <- list()
    for (l in 2:G) { cvec <- numeric(q); cvec[1] <- -1; cvec[l] <- 1
      d <- stats_design_test(rows, rows$unit, m, Z, cvec, tau)
      if (!d$ok) return(NULL)
      tau <- d$tau2; res[[ctx$levels[l]]] <- d }
    res }
  d <- fit(Zg, if (frailty) -1 else 0)
  if (is.null(d)) { out$reason <- "design_not_full_rank"; return(out) }
  out$ok <- TRUE; out$levels_test <- d
  out$test <- list(coef_ref = d[[1]]$theta[1], tau2 = d[[1]]$tau2)
  if (!is.null(Z_extra)) {
    cc <- stats::complete.cases(Z_extra[rows$img + 1L, , drop = FALSE])
    rows_all <- rows; rows <- rows[cc, , drop = FALSE]
    out$covariate_levels <- fit(cbind(Zg[cc, , drop = FALSE], Z_extra[rows$img + 1L, , drop = FALSE]), d[[1]]$tau2)
    rows <- rows_all
  }
  if (availability) {
    i <- rows$img + 1L
    share <- log(pmax(ctx$counts[cbind(i, match(f, ctx$type_labels))], 0.5) / rowSums(ctx$counts)[i])
    if (stats::var(share) > 0) out$availability_levels <- fit(cbind(Zg, share - mean(share)), d[[1]]$tau2)
  }
  out
}

## ---- assembling the results -----------------------------------------------------------------

.cell_table <- function(fits, ctx, availability, covariates) {
  if (length(ctx$levels) > 2L) return(.cell_table_levels(fits, ctx, availability, covariates))
  ok <- vapply(fits, `[[`, logical(1), "ok")
  num <- function(z, nm) if (is.null(z) || !isTRUE(z$ok)) NA_real_ else z[[nm]]
  tab <- do.call(rbind, lapply(fits[ok], function(o) {
    x <- o$test
    row <- data.frame(from = o$from, to = o$to, excess_ref = x$coef_ref, excess_comp = x$coef_comp,
                      excess_difference = x$difference, se = x$se, df = x$df, p_value = x$p, tau2 = x$tau2,
                      stringsAsFactors = FALSE)
    if (availability) {
      row$adjusted_difference <- num(o$availability, "estimate"); row$adjusted_se <- num(o$availability, "se")
      row$adjusted_df <- num(o$availability, "df"); row$adjusted_p_value <- num(o$availability, "p")
    }
    if (!is.null(covariates)) {
      row$covariate_difference <- num(o$covariate, "estimate"); row$covariate_se <- num(o$covariate, "se")
      row$covariate_df <- num(o$covariate, "df"); row$covariate_p_value <- num(o$covariate, "p")
    }
    row
  }))
  if (is.null(tab)) return(NULL)
  tab$p_adj <- stats::p.adjust(tab$p_value, "BH")
  if (availability) tab$adjusted_p_adj <- stats::p.adjust(tab$adjusted_p_value, "BH")
  if (!is.null(covariates)) tab$covariate_p_adj <- stats::p.adjust(tab$covariate_p_value, "BH")
  rownames(tab) <- paste(tab$from, tab$to, sep = "__")
  tab
}

## More than two conditions: one row per pair and level (contrast with the reference level).
.cell_table_levels <- function(fits, ctx, availability, covariates) {
  num <- function(z, l, nm) if (is.null(z) || is.null(z[[l]])) NA_real_ else z[[l]][[nm]]
  tab <- do.call(rbind, lapply(fits[vapply(fits, `[[`, TRUE, "ok")], function(o) {
    do.call(rbind, lapply(names(o$levels_test), function(l) { d <- o$levels_test[[l]]
      row <- data.frame(from = o$from, to = o$to, level = l, excess_ref = d$theta[1],
                        excess_difference = d$estimate, se = d$se, df = d$df, p_value = d$p, tau2 = d$tau2,
                        stringsAsFactors = FALSE)
      if (availability) { row$adjusted_difference <- num(o$availability_levels, l, "estimate")
        row$adjusted_p_value <- num(o$availability_levels, l, "p") }
      if (!is.null(covariates)) { row$covariate_difference <- num(o$covariate_levels, l, "estimate")
        row$covariate_p_value <- num(o$covariate_levels, l, "p") }
      row })) }))
  if (is.null(tab)) return(NULL)
  for (l in unique(tab$level)) { k <- tab$level == l
    tab$p_adj[k] <- stats::p.adjust(tab$p_value[k], "BH")
    if (availability) tab$adjusted_p_adj[k] <- stats::p.adjust(tab$adjusted_p_value[k], "BH") }
  rownames(tab) <- paste(tab$from, tab$to, tab$level, sep = "__")
  tab
}

## Per-image excess (O - E) / n of every pair, for plots and bind(): images x pairs.
.cell_image_excess <- function(fits, ctx) {
  out <- lapply(fits, function(o) {
    v <- rep(NA_real_, ctx$n_images)
    if (!is.null(o$rows) && nrow(o$rows)) v[o$rows$img + 1L] <- (o$rows$O - o$rows$E) / o$rows$n
    v })
  names(out) <- vapply(fits, function(o) paste(o$from, o$to, sep = "__"), "")
  out
}
