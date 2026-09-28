## effect = "excess": the additive (K-scale) effect under random labelling
## (math.tex, Section 18).
##
## For a pair REF -> TARGET in image i, c_j is the number of REF cells among the
## neighbours of candidate cell j (within r, or among its k nearest), candidates
## being the non-REF cells (all cells for a self-pair), and
##   O_i = sum over TARGET cells of c_j,   L_i = sum over candidates of c_j,
##   Q_i = sum over candidates of c_j^2,   p_i = n_TARGET / M_i.
## Under random labelling of the candidates E(O_i) = p_i L_i and
##   v_i = Var(O_i) = p_i (1 - p_i) M_i / (M_i - 1) (Q_i - L_i^2 / M_i).
## The group parameter delta_g is the mean number of extra REF neighbours per
## TARGET cell beyond random labelling, lambda_REF (K_obs - K_RL) for the radius
## graph: O_i = p_i L_i + delta_g n_i + e_i. Recruiting a fraction f of TARGET
## cells next to REF cells gives delta = f whatever the density or abundance.
##
## Working covariance within a unit u: diag(v_i) + tau2 n n^T, the rank-one
## structure of the frailty GEE with I_i = n_i (Section 17), so the unit score is
##   r_u = sum_i n_i (O_i - p_i L_i - delta n_i) / v_i,  J_u = sum_i n_i^2 / v_i,
## the weight a_u = 1 / (1 + tau2 J_u), tau2 is Paule-Mandel and the variance the
## closed-form CR2 (frailty_cr2) with Satterthwaite df. The model is linear in
## delta, so the estimate is explicit.

excess_image_data <- function(ctx, from_code, to_code) {
  tot <- ctx$totals; sq <- ctx$out_sq_totals; counts <- ctx$counts
  from <- from_code + 1L; to <- to_code + 1L; self <- from == to
  knn <- ctx$family == "binomial"
  nA <- counts[, from]; nB <- counts[, to]; N <- rowSums(counts)
  cand <- if (self) rep(TRUE, ncol(counts)) else seq_len(ncol(counts)) != from
  # tot[a, b, i] = sum over the image's b cells of their a neighbours; the radius
  # totals count each cell as its own neighbour, which matters for a self-pair only
  O <- as.numeric(tot[from, to, ]); L <- as.numeric(apply(tot[from, cand, , drop = FALSE], 3, sum))
  if (self && !knn) { O <- O - nA; L <- L - nA }
  Q <- as.numeric(apply(sq[cand, from, , drop = FALSE], 3, sum))
  M <- if (self) N else N - nA
  p <- if (self) pmax(nA - 1, 0) / pmax(M - 1, 1) else nB / pmax(M, 1)
  v <- p * (1 - p) * M / pmax(M - 1, 1) * pmax(Q - L^2 / pmax(M, 1), 0)
  if (!is.null(ctx$phi_inflation)) v <- v * ctx$phi_inflation[to, ]
  im <- data.frame(img = seq_len(nrow(counts)) - 1L, O = O, E = p * L, n = nB, v = v)
  im[im$n > 0 & im$v > 0, , drop = FALSE]
}

excess_fit_group <- function(d, unit, tau2) {
  J <- as.numeric(tapply(d$n^2 / d$v, unit, sum)); a <- 1 / (1 + tau2 * J)
  s <- as.numeric(tapply(d$n * (d$O - d$E) / d$v, unit, sum))
  delta <- sum(a * s) / sum(a * J)
  r <- as.numeric(tapply(d$n * (d$O - d$E - delta * d$n) / d$v, unit, sum))
  list(beta = delta, r = r, J = J, a = a, B = sum(a * J))
}

excess_pearson <- function(d, unit, group, tau2) {
  X <- 0
  for (g in 0:1) { k <- group == g
    f <- excess_fit_group(d[k, , drop = FALSE], droplevels(unit[k]), tau2); X <- X + sum(f$a * f$r^2 / f$J) }
  X
}

# Paule-Mandel, as frailty_tau2; tau2 is on the scale of delta^2 (neighbours per cell)
excess_tau2 <- function(d, unit, group) {
  df <- nlevels(unit) - 2
  if (df < 1 || excess_pearson(d, unit, group, 0) <= df) return(0)
  up <- 1e-4
  while (excess_pearson(d, unit, group, up) > df && up < 1e6) up <- up * 4
  stats::uniroot(function(t) excess_pearson(d, unit, group, t) - df, c(0, up), tol = 1e-10)$root
}

excess_outcome <- function(ctx, f, t) {
  im <- excess_image_data(ctx, ctx$type_index[[f]], ctx$type_index[[t]])
  unit <- factor(ctx$image_cluster[im$img + 1L]); group <- ctx$image_group[im$img + 1L]
  for (g in 0:1) if (length(unique(unit[group == g])) < 2L)
    return(skip_pair(f, t, "one_patient_per_group", paste0(
      "Skipping pair ", f, "__", t, ": condition '", ctx$levels[g + 1L],
      "' has fewer than two clusters with both cell types; CR2 needs at least two per group.")))
  tau2 <- if (ctx$frailty) excess_tau2(im, unit, group) else 0
  fits <- lapply(0:1, function(g) { k <- group == g
    excess_fit_group(im[k, , drop = FALSE], droplevels(unit[k]), tau2) })
  cr <- lapply(fits, frailty_cr2)
  v_hat <- cr[[1]]$V + cr[[2]]$V; df <- (cr[[1]]$EV + cr[[2]]$EV)^2 / (cr[[1]]$trsq + cr[[2]]$trsq)
  diff <- fits[[2]]$beta - fits[[1]]$beta
  list(from = f, to = t, condition_ref = ctx$levels[1], condition_comp = ctx$levels[2],
       coef_ref = fits[[1]]$beta, coef_comp = fits[[2]]$beta, excess_difference = diff,
       p_value = 2 * stats::pt(-abs(diff / sqrt(v_hat)), df),
       estimator = "gee", family = ctx$family, mle_would_skip = FALSE, mle_skip_reason = NA_character_,
       .frailty = c(tau2 = tau2, df = df, se = sqrt(v_hat)),
       .images = if (isTRUE(ctx$moderate)) list(im = im, unit = unit, group = group))
}

# moderate = TRUE with effect = "excess": as frailty_moderate(). Each pair keeps its own
# estimate and CR2 variance V (Satterthwaite df d); the prior for V / Vm, with Vm the working-model
# variance sum_g 1 / B_g at the median tau2 across pairs, is Smyth's scaled inverse chi-square
# (fit_f_dist, no trend), and the test uses the posterior variance on d + d0 df.
excess_moderate <- function(outcomes, ctx) {
  ok <- which(vapply(outcomes, function(o) is.null(o$reason), logical(1)))
  if (length(ok) < 3L) stop("moderate = TRUE needs at least three pairs that can be fitted.", call. = FALSE)
  tau2 <- stats::median(vapply(outcomes[ok], function(o) o$.frailty[["tau2"]], numeric(1)))
  R <- t(vapply(outcomes[ok], function(o) { x <- o$.images
    B <- vapply(0:1, function(g) { k <- x$group == g
      excess_fit_group(x$im[k, , drop = FALSE], droplevels(x$unit[k]), tau2)$B }, numeric(1))
    c(V = o$.frailty[["se"]]^2, d = o$.frailty[["df"]], Vm = sum(1 / B)) }, numeric(3)))
  prior <- fit_f_dist(R[, "V"] / R[, "Vm"], R[, "d"], rep(0, nrow(R)))
  s0 <- prior$s20[1]; d0 <- prior$d0
  ratio <- R[, "V"] / R[, "Vm"]
  post <- if (is.finite(d0)) (d0 * s0 + R[, "d"] * ratio) / (d0 + R[, "d"]) else rep(s0, nrow(R))
  df <- pmin(R[, "d"] + d0, sum(R[, "d"]))
  for (j in seq_along(ok)) { i <- ok[j]; o <- outcomes[[i]]
    se <- unname(sqrt(post[j] * R[j, "Vm"]))
    o$p_value <- unname(2 * stats::pt(-abs(o$excess_difference / se), df[j]))
    o$.frailty <- c(tau2 = o$.frailty[["tau2"]], df = unname(df[j]), se = se)
    outcomes[[i]] <- o }
  attr(outcomes, "frailty_prior") <- c(tau2 = unname(tau2), d0 = unname(d0), s0 = unname(s0))
  outcomes
}
