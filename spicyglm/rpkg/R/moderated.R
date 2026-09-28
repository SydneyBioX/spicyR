## test = "moderated": image-level summaries compared between conditions with an
## empirical-Bayes moderated t-test across pairs (limma's eBayes(trend = TRUE)).

# Per-image summary of one pair: theta = log((O + 1/2) / E) for the Poisson
# designs, the log odds ratio logit((O + 1/2) / (T + 1)) - logit(E / T) for the
# binomial ones; NA where the image has no reference cell for the pair.
image_summary <- function(md, n_images, family, k) {
  im <- factor(md$image, levels = seq_len(n_images) - 1L)
  O <- as.numeric(tapply(as.numeric(md$n), im, sum))
  if (family == "binomial") {
    trials <- if (!is.null(md$trials)) as.numeric(md$trials) else rep(k, length(md$n))
    Tr <- as.numeric(tapply(trials, im, sum))
    E <- as.numeric(tapply(trials * md$p0, im, sum))
    theta <- stats::qlogis((O + 0.5) / (Tr + 1)) - stats::qlogis(E / Tr)
  } else {
    E <- as.numeric(tapply(md$density, im, sum))
    theta <- log((O + 0.5) / E)
  }
  theta[!is.finite(theta)] <- NA_real_
  list(O = O, theta = theta)
}

# image_summary() of a plain (homogeneous) design from the per-image totals,
# over the images that poisson_model_data / binomial_model_data would keep.
totals_summary <- function(totals, counts, ctx, f, t) {
  O <- totals[t, f, ]
  n_from <- as.numeric(counts[, f]); n_to <- as.numeric(counts[, t])
  if (ctx$family == "binomial") {
    n_all <- rowSums(counts)
    keep <- n_from > 0 & n_to > 0 & n_all > ctx$k & n_to < n_all
    Tr <- n_from * ctx$k
    theta <- stats::qlogis((O + 0.5) / (Tr + 1)) - stats::qlogis(n_to / n_all)
  } else {
    keep <- n_from > 0 & n_to > 0
    theta <- log((O + 0.5) / (n_from * n_to / ctx$areas * pi * ctx$r^2))
  }
  O[!keep] <- NA_real_
  theta[!keep | !is.finite(theta)] <- NA_real_
  list(O = O, theta = theta)
}

# Slopes of theta on the log abundance of the two types, shared by all pairs and
# estimated within pair x condition (a fixed effect for each), so that neither a
# pair's level nor its condition effect enters them: an ANCOVA with a common slope.
# Only images with at least 5 observed neighbours are used.
density_slopes <- function(summaries, pairs, counts, image_group) {
  rows <- lapply(seq_along(pairs), function(i) {
    s <- summaries[[i]]; if (is.null(s)) return(NULL)
    ok <- !is.na(s$theta) & s$O >= 5
    if (sum(ok) < 3L) return(NULL)
    X <- cbind(log(pmax(counts[ok, pairs[[i]][1]], 1)), log(pmax(counts[ok, pairs[[i]][2]], 1)), s$theta[ok])
    g <- image_group[ok]
    for (k in unique(g)) X[g == k, ] <- sweep(X[g == k, , drop = FALSE], 2L, colMeans(X[g == k, , drop = FALSE]))
    X
  })
  X <- do.call(rbind, rows)
  if (is.null(X) || nrow(X) < 10L) return(c(0, 0))
  b <- qr.coef(qr(X[, 1:2]), X[, 3])
  b[is.na(b)] <- 0                            # collinear only if every pair is a self pair
  b
}

# limma's trigammaInverse: solve trigamma(y) = x by Newton's method.
trigamma_inverse <- function(x) {
  y <- 0.5 + 1 / x
  for (iter in 1:50) {
    tri <- trigamma(y)
    dif <- tri * (1 - tri / x) / psigamma(y, deriv = 2L)
    y <- y + dif
    if (max(-dif / y) < 1e-8) break
  }
  y
}

# limma's fitFDist with a trend on `covariate`: prior df d0 and prior variance s0^2.
fit_f_dist <- function(s2, d, covariate) {
  m <- stats::median(s2); if (m <= 0) m <- 1
  s2 <- pmax(s2, 1e-5 * m)
  e <- log(s2) - digamma(d / 2) + log(d / 2)
  spline_df <- min(1L + (length(e) >= 3L) + (length(e) >= 6L) + (length(e) >= 30L),
                   length(unique(covariate)))
  if (spline_df < 2L) {
    fitted <- rep(mean(e), length(e)); evar <- sum((e - mean(e))^2) / (length(e) - 1)
  } else {
    design <- splines::ns(covariate, df = spline_df, intercept = TRUE)
    fit <- stats::lm.fit(design, e)
    fitted <- fit$fitted.values
    evar <- if (fit$df.residual > 0) sum(fit$residuals^2) / fit$df.residual else 0
  }
  evar <- evar - mean(trigamma(d / 2))
  if (evar > 0) {
    d0 <- 2 * trigamma_inverse(evar)
    s20 <- exp(fitted + digamma(d0 / 2) - log(d0 / 2))
  } else {
    d0 <- Inf; s20 <- if (spline_df < 2L) rep(mean(s2), length(s2)) else exp(fitted)
  }
  list(d0 = d0, s20 = s20)
}

moderated_tests <- function(ctx, pairs, density_adjust, counts) {
  n_images <- length(ctx$image_group)
  eff <- EFFECT_COLUMNS[[ctx$family]]

  outcomes <- vector("list", length(pairs))
  summaries <- vector("list", length(pairs))
  # the plain designs' summaries come from one pass over the neighbour graph
  totals <- if (is.null(ctx$parent) && is.null(ctx$sigma))
    dataset_pair_neighbour_totals(ctx$data, ctx$family == "binomial", ncol(counts))
  for (i in seq_along(pairs)) {
    f <- pairs[[i]][1]; t <- pairs[[i]][2]
    if (!is.null(totals)) {
      s <- totals_summary(totals, counts, ctx, ctx$type_index[[f]] + 1L, ctx$type_index[[t]] + 1L)
      images <- which(!is.na(s$theta)) - 1L
    } else {
      md <- pair_model_data(ctx, f, t)
      images <- md$image
    }
    present <- unique(ctx$image_group[images + 1L])
    if (length(present) < 2L) {
      d <- diagnose_missing(f, t, setdiff(0:1, present), ctx)
      outcomes[[i]] <- skip_pair(f, t, d$reason, d$message)
      next
    }
    summaries[[i]] <- if (!is.null(totals)) s else image_summary(md, n_images, ctx$family, ctx$k)
  }
  codes <- lapply(pairs, function(p) c(ctx$type_index[[p[1]]], ctx$type_index[[p[2]]]) + 1L)
  mod <- moderate(summaries, codes, counts, ctx$image_group, ctx$image_cluster, density_adjust)

  for (i in mod$too_few) {
    f <- pairs[[i]][1]; t <- pairs[[i]][2]
    outcomes[[i]] <- skip_pair(f, t, "one_patient_per_group", paste0(
      "Skipping pair ", f, "__", t, ": a condition has fewer than two ",
      "subjects with data for this pair; the moderated test needs at least two per group."))
  }
  tab <- mod$table
  for (j in seq_len(nrow(tab))) {
    i <- tab$index[j]
    out <- list(from = pairs[[i]][1], to = pairs[[i]][2],
                condition_ref = ctx$levels[1], condition_comp = ctx$levels[2],
                coef_ref = tab$mean_ref[j], coef_comp = tab$mean_comp[j],
                log_effect = tab$estimate[j], effect = exp(tab$estimate[j]),
                p_value = tab$p_value[j], estimator = "moderated", family = ctx$family,
                mle_would_skip = FALSE, mle_skip_reason = NA_character_)
    names(out)[names(out) == "log_effect"] <- eff[1]
    names(out)[names(out) == "effect"] <- eff[2]
    outcomes[[i]] <- out
  }
  tab$from <- vapply(pairs[tab$index], `[`, "", 1L)
  tab$to <- vapply(pairs[tab$index], `[`, "", 2L)
  attr(outcomes, "moderation") <- list(
    prior_df = mod$prior_df,
    density_slopes = if (density_adjust) stats::setNames(mod$slopes, c("from", "to")),
    pairs = tab[, c("from", "to", "mean_log_count", "s2", "df", "prior_s2", "posterior_s2",
                    "df_total", "t")])
  outcomes
}

# The moderated test on per-image summaries (NULL for pairs not fitted). Units
# are the clusters in `image_unit`; each unit's summary is the mean over its
# images, and a condition is a unit-level label.
moderate <- function(summaries, codes, counts, image_group, image_unit, density_adjust) {
  unit_group <- as.integer(tapply(image_group, image_unit, `[`, 1L))
  slopes <- if (density_adjust) density_slopes(summaries, codes, counts, image_group) else c(0, 0)

  idx <- which(!vapply(summaries, is.null, logical(1)))
  U <- matrix(NA_real_, length(idx), length(unit_group)); A <- numeric(length(idx))
  for (j in seq_along(idx)) {
    s <- summaries[[idx[j]]]; theta <- s$theta
    if (density_adjust)
      theta <- theta - slopes[1] * log(pmax(counts[, codes[[idx[j]]][1]], 1)) -
        slopes[2] * log(pmax(counts[, codes[[idx[j]]][2]], 1))
    u <- as.numeric(tapply(theta, image_unit, mean, na.rm = TRUE))
    U[j, ] <- ifelse(is.finite(u), u, NA_real_)
    A[j] <- mean(log1p(s$O[!is.na(s$theta)]))
  }
  in0 <- !is.na(U) & rep(unit_group == 0L, each = nrow(U))
  in1 <- !is.na(U) & rep(unit_group == 1L, each = nrow(U))
  n0 <- rowSums(in0); n1 <- rowSums(in1)
  ok <- n0 >= 2 & n1 >= 2
  if (sum(ok) < 3L)
    stop("test = 'moderated' needs at least three pairs that can be fitted; use test = 'glm'.",
         call. = FALSE)
  U0 <- ifelse(in0, U, 0); U1 <- ifelse(in1, U, 0)
  m0 <- rowSums(U0) / n0; m1 <- rowSums(U1) / n1
  rss <- rowSums(ifelse(in0, (U - m0)^2, 0)) + rowSums(ifelse(in1, (U - m1)^2, 0))
  d <- n0 + n1 - 2; s2 <- rss / d; cc <- 1 / n0 + 1 / n1

  prior <- fit_f_dist(s2[ok], d[ok], A[ok])
  s2_post <- if (is.finite(prior$d0)) (prior$d0 * prior$s20 + d[ok] * s2[ok]) / (prior$d0 + d[ok])
             else prior$s20
  df_total <- pmin(prior$d0 + d[ok], sum(d[ok]))
  est <- (m1 - m0)[ok]
  t_stat <- est / sqrt(s2_post * cc[ok])
  list(table = data.frame(index = idx[ok], mean_ref = m0[ok], mean_comp = m1[ok], estimate = est,
                          t = t_stat, p_value = 2 * stats::pt(-abs(t_stat), df_total),
                          mean_log_count = A[ok], s2 = s2[ok], df = d[ok], prior_s2 = prior$s20,
                          posterior_s2 = s2_post, df_total = df_total),
       too_few = idx[!ok], prior_df = prior$d0, slopes = slopes)
}
