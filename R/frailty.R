## frailty = TRUE: the cell-level GLM with an image (subject) frailty and the
## random-labelling within-image variance (math.tex, Section 17).
##
## For a unit u (a subject; its images i) the working covariance of the image
## totals is diag(phi_i I_i) + tau2 I I^T, with I_i = dm_i / dbeta and phi_i the
## random-labelling variance of O_i over its working variance. The GEE score of a
## group is sum_u a_u r_u with
##   r_u = sum_i (O_i - m_i) / phi_i,  J_u = sum_i I_i / phi_i,  a_u = 1 / (1 + tau2 J_u),
## i.e. the cell-level GLM with prior weight a_u / phi_i on the cells of image i.
## The working covariance is diagonal plus rank one within a unit, so CR2 has
## the closed form A_u = (1 - h_u)^{-1/2}, h_u = a_u J_u / B.

# phi_i for every image of a pair, from the per-image neighbour totals and their
# sums of squares ([to, from, image] arrays), under random labelling of the
# candidates: the non-REF cells, or all other cells for a self-pair.
frailty_phi <- function(tot, sq, counts, from, to, knn) {
  nA <- counts[, from]; nB <- counts[, to]; N <- rowSums(counts)
  self <- from == to
  cand <- if (self) rep(TRUE, dim(tot)[1]) else seq_len(dim(tot)[1]) != from
  L <- as.numeric(colSums(tot[cand, from, , drop = FALSE]))
  if (self && !knn) L <- L - nA          # radius totals count each cell as its own neighbour
  Q <- as.numeric(colSums(sq[cand, from, , drop = FALSE]))
  Np <- if (self) N - 1 else N - nA
  p <- if (self) (nA - 1) / Np else nB / Np
  v <- p * (1 - p) * Np / pmax(Np - 1, 1) * pmax(Q - L^2 / Np, 0)
  phi <- if (knn) v / (L * p * (1 - p)) else v / (p * L)
  phi[!is.finite(phi) | phi <= 0] <- 1
  phi
}

# phi for the weighted designs. Kontextual: random labelling of the context cells
# (the context and its local densities do not change under relabelling).
# Inhomogeneous: the TARGET weights depend on its own estimated intensity, so the
# null is an inhomogeneous Poisson process for the TARGET; by Campbell's theorem
# Var(O) = int c^2 lambda, estimated without bias by sum over TARGET cells of c_b^2.
frailty_phi_weighted <- function(ctx, from_code, to_code) {
  design <- if (!is.null(ctx$sigma)) 2L else if (ctx$family == "binomial") 1L else 0L
  S <- dataset_weighted_phi_sums(ctx$data, from_code, to_code, design, isTRUE(ctx$edge_correct))
  s1 <- S[1, ]; s2 <- S[2, ]; Np <- S[3, ]
  if (design == 2L) { phi <- s2 / s1 }
  else {
    nB <- ctx$counts[, to_code + 1L]; self <- from_code == to_code
    p <- if (self) (nB - 1) / pmax(Np - 1, 1) else nB / Np
    v <- p * (1 - p) * Np / pmax(Np - 1, 1) * pmax(s2 - s1^2 / Np, 0)
    phi <- if (design == 1L) v / (s1 * p * (1 - p)) else v / (p * s1)
  }
  phi[!is.finite(phi) | phi <= 0] <- 1
  phi
}

frailty_phi_all <- function(ctx, from_code, to_code) {
  if (!is.null(ctx$parent) || !is.null(ctx$sigma)) frailty_phi_weighted(ctx, from_code, to_code)
  else frailty_phi(ctx$totals, ctx$sq_totals, ctx$counts, from_code + 1L, to_code + 1L, ctx$family == "binomial")
}

# Label clustering (TARGET cells clustered among themselves) inflates Var(O) above
# random labelling. The spatial HAC estimate captures it but is noisy image by image,
# so it is pooled: per (TARGET type, image), the median over REF types of
# HAC / E[HAC under random labelling], floored at 1. The HAC is compared with its own
# null expectation, not with phi: a bandwidth-limited HAC keeps only the nearby part of
# the negative finite-population covariance that phi subtracts in full, so HAC / phi
# exceeds 1 under random labelling whenever the REF type is common.
# Returns an n_types x n_images matrix.
frailty_label_clustering <- function(ctx, pairs) {
  n_types <- ncol(ctx$counts); n_img <- nrow(ctx$counts)
  design <- if (!is.null(ctx$parent)) (if (ctx$family == "binomial") 1L else 0L) else if (ctx$family == "binomial") 4L else 3L
  h <- ctx$hac_h
  # one neighbour pass per REF type for all its non-self TARGETs (dataset_hac_phi_sums_ref);
  # self-pairs, whose candidates include the REF type, use the per-pair sums
  codes <- t(vapply(pairs, function(p) c(ctx$type_index[[p[1]]], ctx$type_index[[p[2]]]), integer(2)))
  by_ref <- lapply(split(seq_len(nrow(codes)), codes[, 1]), function(k) {
    f <- codes[k[1], 1]
    if (any(codes[k, 2] != f)) dataset_hac_phi_sums_ref(ctx$data, f, design, h, n_types) })
  ratios <- lapply(seq_len(nrow(codes)), function(i) { f <- codes[i, 1]; t <- codes[i, 2]
    if (f == t) { H <- dataset_hac_phi_sums(ctx$data, f, t, design, h); V <- H[1, ]; Np <- H[3, ]; G <- H[4, ] }
    else { H <- by_ref[[as.character(f)]]; V <- H[t + 1L, ]; Np <- H[n_types + 2L, ]; G <- H[n_types + 3L, ] }
    nB <- ctx$counts[, t + 1L]; pr <- if (f == t) (nB - 1) / pmax(Np - 1, 1) else nB / Np
    ratio <- V / (pr * (1 - pr) * G)
    O <- dataset_pair_image_totals(ctx, f, t); ratio[O < 5] <- NA       # sparse images: too noisy to pool
    list(to = t + 1L, ratio = ratio) })
  out <- matrix(1, n_types, n_img)
  for (t in unique(vapply(ratios, `[[`, integer(1), "to"))) {
    R <- do.call(rbind, lapply(Filter(function(x) x$to == t, ratios), `[[`, "ratio"))
    R[!is.finite(R) | R <= 0] <- NA
    k <- apply(R, 2, function(z) if (all(is.na(z))) 1 else stats::median(z, na.rm = TRUE))
    out[t, ] <- pmax(k, 1)
  }
  out
}

# O_i of a pair for every image, from the one-pass neighbour totals (radius or k-NN).
dataset_pair_image_totals <- function(ctx, f, t) {
  if (is.null(ctx$pair_totals)) return(rep(Inf, nrow(ctx$counts)))
  as.numeric(ctx$pair_totals[t + 1L, f + 1L, ])
}

frailty_mean <- function(d, beta, binomial) {
  if (!binomial) { m <- d$E * exp(beta); return(list(m = m, I = m, dI = m)) }
  pr <- stats::plogis(beta + stats::qlogis(d$p0)); m <- d$Tr * pr; I <- m * (1 - pr)
  list(m = m, I = I, dI = I * (1 - 2 * pr))
}

# One group's Firth-adjusted GEE fit at fixed tau2: beta and the unit pieces.
frailty_fit_group <- function(d, unit, tau2, binomial) {
  beta <- if (!binomial) log((sum(d$O) + 0.5) / sum(d$E))
          else stats::qlogis((sum(d$O) + 0.5) / (sum(d$Tr) + 1)) - stats::qlogis(sum(d$Tr * d$p0) / sum(d$Tr))
  pieces <- function(beta) { mf <- frailty_mean(d, beta, binomial)
    r <- as.numeric(tapply((d$O - mf$m) / d$phi, unit, sum)); J <- as.numeric(tapply(mf$I / d$phi, unit, sum))
    dJ <- as.numeric(tapply(mf$dI / d$phi, unit, sum)); a <- 1 / (1 + tau2 * J)
    list(r = r, J = J, dJ = dJ, a = a, B = sum(a * J)) }
  for (it in 1:100) {
    p <- pieces(beta)
    step <- (sum(p$a * p$r) + 0.5 * sum(p$a * p$dJ) / p$B) / p$B
    if (!is.finite(step)) break
    step <- max(min(step, 2), -2); beta <- beta + step
    if (abs(step) < 1e-10) break
  }
  c(list(beta = beta), pieces(beta))
}

# CR2 variance of a group's beta and the two traces of its Satterthwaite df.
frailty_cr2 <- function(f) {
  n <- length(f$r); h <- f$a * f$J / f$B; cc <- f$a^2 / f$B^2
  H <- outer(f$J, f$a) / f$B
  E <- diag(1 / sqrt(1 - h), n) %*% (diag(n) - H)
  Sig <- E %*% (diag(f$J / f$a, n) %*% t(E))
  MS <- cc * Sig
  list(V = sum(cc * f$r^2 / (1 - h)), EV = sum(diag(MS)), trsq = sum(MS * t(MS)))
}

frailty_pearson <- function(d, unit, group, tau2, binomial) {
  X <- 0
  for (g in 0:1) { k <- group == g
    f <- frailty_fit_group(d[k, , drop = FALSE], droplevels(unit[k]), tau2, binomial); X <- X + sum(f$a * f$r^2 / f$J) }
  X
}

# Paule-Mandel: the tau2 at which the Pearson statistic equals its degrees of freedom.
frailty_tau2 <- function(d, unit, group, binomial) {
  df <- nlevels(unit) - 2
  if (df < 1 || frailty_pearson(d, unit, group, 0, binomial) <= df) return(0)
  up <- 1
  while (frailty_pearson(d, unit, group, up, binomial) > df && up < 1e4) up <- up * 4
  stats::uniroot(function(t) frailty_pearson(d, unit, group, t, binomial) - df, c(0, up), tol = 1e-6)$root
}

# The frailty test for one pair. `images` is a data frame with one row per image
# holding img (0-based), O, and E (Poisson) or Tr and p0 (binomial), plus phi.
frailty_pair <- function(images, ctx) {
  binomial <- ctx$family == "binomial"
  unit <- factor(ctx$image_cluster[images$img + 1L]); group <- ctx$image_group[images$img + 1L]
  tau2 <- frailty_tau2(images, unit, group, binomial)
  fits <- lapply(0:1, function(g) { k <- group == g
    frailty_fit_group(images[k, , drop = FALSE], droplevels(unit[k]), tau2, binomial) })
  cr <- lapply(fits, frailty_cr2)
  list(beta = c(fits[[1]]$beta, fits[[2]]$beta), v_hat = cr[[1]]$V + cr[[2]]$V,
       df = (cr[[1]]$EV + cr[[2]]$EV)^2 / (cr[[1]]$trsq + cr[[2]]$trsq), tau2 = tau2)
}

# Collapse a pair's cell-level model data to the image rows frailty_pair() needs.
frailty_images <- function(md, ctx, from_code, to_code) {
  keep <- if (!is.null(md$density)) md$density > 0 else rep(TRUE, length(md$n))
  im <- factor(md$image[keep], levels = sort(unique(md$image[keep])))
  img <- as.integer(levels(im))
  O <- as.numeric(tapply(md$n[keep], im, sum))
  phi <- frailty_phi_all(ctx, from_code, to_code)
  if (!is.null(ctx$phi_inflation)) phi <- phi * ctx$phi_inflation[to_code + 1L, ]
  phi <- phi[img + 1L]
  if (ctx$family == "binomial") {
    trials <- if (!is.null(md$trials)) as.numeric(md$trials[keep]) else rep(ctx$k, sum(keep))
    Tr <- as.numeric(tapply(trials, im, sum)); p0 <- as.numeric(tapply(md$p0[keep], im, `[`, 1))
    data.frame(img = img, O = O, Tr = Tr, p0 = p0, phi = phi)
  } else {
    data.frame(img = img, O = O, E = as.numeric(tapply(md$density[keep], im, sum)), phi = phi)
  }
}

frailty_outcome <- function(md, ctx, f, t, b) {
  from_code <- ctx$type_index[[f]]; to_code <- ctx$type_index[[t]]
  images <- frailty_images(md, ctx, from_code, to_code)
  unit_group <- tapply(ctx$image_group[images$img + 1L], ctx$image_cluster[images$img + 1L], `[`, 1)
  for (g in 0:1) if (sum(unit_group == g) < 2L)
    return(skip_pair(f, t, "one_patient_per_group", paste0(
      "Skipping pair ", f, "__", t, ": condition '", ctx$levels[g + 1L],
      "' has fewer than two clusters with data for this pair; CR2 needs at least two per group.")))
  fit <- frailty_pair(images, ctx)
  log_effect <- fit$beta[2] - fit$beta[1]
  eff <- EFFECT_COLUMNS[[ctx$family]]
  out <- list(from = f, to = t, condition_ref = ctx$levels[1], condition_comp = ctx$levels[2],
              coef_ref = fit$beta[1], coef_comp = fit$beta[2],
              log_effect = log_effect, effect = exp(log_effect),
              p_value = 2 * stats::pt(-abs(log_effect / sqrt(fit$v_hat)), fit$df),
              estimator = "firth", family = ctx$family,
              mle_would_skip = !is.null(b$reason),
              mle_skip_reason = if (is.null(b$reason)) NA_character_ else b$reason,
              .frailty = c(tau2 = fit$tau2, df = fit$df, se = sqrt(fit$v_hat)),
              .images = if (isTRUE(ctx$moderate)) images)
  names(out)[names(out) == "log_effect"] <- eff[1]
  names(out)[names(out) == "effect"] <- eff[2]
  out
}

# moderate = TRUE: every pair refitted at the median tau2; the CR2 variance V
# (Satterthwaite df d) is shrunk toward the frailty-model variance Vm with the
# prior V / Vm ~ s0 * (scaled inverse chi-square on d0 df), fitted across pairs.
frailty_moderate <- function(outcomes, ctx) {
  ok <- which(vapply(outcomes, function(o) is.null(o$reason), logical(1)))
  if (length(ok) < 3L) stop("moderate = TRUE needs at least three pairs that can be fitted.", call. = FALSE)
  tau2 <- stats::median(vapply(outcomes[ok], function(o) o$.frailty[["tau2"]], numeric(1)))
  binomial <- ctx$family == "binomial"
  eff <- EFFECT_COLUMNS[[ctx$family]]
  # each pair keeps its own estimate and CR2 (weights at its own tau2); the prior
  # variance is its frailty-model variance at the common tau2
  R <- t(vapply(outcomes[ok], function(o) { im <- o$.images
    unit <- factor(ctx$image_cluster[im$img + 1L]); group <- ctx$image_group[im$img + 1L]
    fits <- lapply(0:1, function(g) { k <- group == g
      frailty_fit_group(im[k, , drop = FALSE], droplevels(unit[k]), tau2, binomial) })
    c(V = o$.frailty[["se"]]^2, d = o$.frailty[["df"]], Vm = 1 / fits[[1]]$B + 1 / fits[[2]]$B) }, numeric(3)))
  prior <- fit_f_dist(R[, "V"] / R[, "Vm"], R[, "d"], rep(0, nrow(R)))   # constant covariate: no trend
  s0 <- prior$s20[1]; d0 <- prior$d0
  ratio <- R[, "V"] / R[, "Vm"]
  post <- if (is.finite(d0)) (d0 * s0 + R[, "d"] * ratio) / (d0 + R[, "d"]) else rep(s0, nrow(R))
  df <- pmin(R[, "d"] + d0, sum(R[, "d"]))
  for (j in seq_along(ok)) { i <- ok[j]; o <- outcomes[[i]]
    se <- unname(sqrt(post[j] * R[j, "Vm"]))
    o$p_value <- unname(2 * stats::pt(-abs(o[[eff[1]]] / se), df[j]))
    o$.frailty <- c(tau2 = o$.frailty[["tau2"]], df = unname(df[j]), se = se)
    outcomes[[i]] <- o }
  attr(outcomes, "frailty_prior") <- c(tau2 = unname(tau2), d0 = unname(d0), s0 = unname(s0))
  outcomes
}
