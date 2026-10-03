## Plain-R reference for spicyR Cell, from ProjectSpicyRCell R/reference/excess_reference.R (written from the
## Supplementary's general definitions, not the package's closed forms). Pairs are in the core's (counted, centre)
## order: ref_* (from = counted, to = centre).
ref_neighbours <- function(x, y, r) { D <- as.matrix(stats::dist(cbind(x, y))); A <- D <= r; diag(A) <- FALSE; list(D = D, adj = A) }

## For pair A -> B in one image: candidates, c_k, y_k and the moments of Proposition 1 (self-pairs: Section 2.3).
ref_image_counts <- function(type, nb, A, B) {
  self <- A == B
  cand <- if (self) rep(TRUE, length(type)) else type != A
  ck <- rowSums(nb$adj[cand, type == A, drop = FALSE])          # A cells among the neighbours of each candidate
  yk <- as.numeric(type[cand] == B)
  M <- sum(cand); nA <- sum(type == A); N <- length(type)
  n <- sum(yk)
  p <- if (self) (nA - 1) / (N - 1) else n / M
  O <- sum(ck * yk); L <- sum(ck); L2 <- sum(ck^2)
  sigma2 <- p * (1 - p) * M / (M - 1)
  list(self = self, cand = cand, ck = ck, yk = yk, M = M, n = n, p = p, O = O, L = L, L2 = L2, E = p * L,
       sigma2 = sigma2, v = sigma2 * (L2 - L^2 / M))
}

## The allocation design: every candidate's score is g_k = 1{c_k >= 1} (any A cell among its neighbours), and
## n = n_B - E, so (O - E) / n is the extra fraction of B cells next to A. A self-pair's E is the exact
## random-labelling expectation (n_A / N) sum_b [1 - C(N - 1 - d_b, n_A - 1) / C(N - 1, n_A - 1)].
ref_image_counts_alloc <- function(type, nb, A, B) {
  ic <- ref_image_counts(type, nb, A, B)
  ic$ck <- as.numeric(ic$ck >= 1)
  ic$O <- sum(ic$ck * ic$yk); ic$L <- ic$L2 <- sum(ic$ck)
  ic$E <- if (ic$self) {
    N <- length(type); nA <- sum(type == A); d <- rowSums(nb$adj)
    if (nA < 2) 0 else nA / N * sum(1 - exp(lchoose(N - 1 - d, nA - 1) - lchoose(N - 1, nA - 1)))
  } else ic$p * ic$L
  ic$v <- ic$sigma2 * (ic$L2 - ic$L^2 / ic$M)
  ic$n_to <- ic$n; ic$n <- ic$n - ic$E
  ic
}

## Exact random-labelling mean and variance of O by enumerating every labelling (tiny images only): checks
## Proposition 1 itself.
ref_rl_enumerate <- function(ck, n) {
  S <- utils::combn(length(ck), n); O <- apply(S, 2, function(s) sum(ck[s]))
  c(mean = mean(O), var = mean(O^2) - mean(O)^2)
}

## ---- Section 3.3: the label-clustering factor -------------------------------------------------------------------

## Per-image factor psi~(A, B) = HAC / (sigma^2 R), with the explicit projector onto span(1, c).
## Returns NA when the factor does not enter the pooling (R <= 0, HAC <= 0, not finite, fewer than five
## (A, B) neighbour pairs, fewer than three candidates).
## self_count = "package" reproduces a known deviation of the package (found 2 Oct 2026): for a self-pair its
## five-pair check counts each A cell as its own neighbour (O + n_A), because it reads the raw radius totals.
ref_psi_tilde <- function(ic, nb, omega, self_count = c("spec", "package")) {
  O5 <- if (match.arg(self_count) == "package" && ic$self) ic$O + ic$n else ic$O
  if (ic$M < 3 || O5 < 5) return(NA_real_)
  X <- cbind(1, ic$ck); qx <- qr(X); Qx <- qr.Q(qx)[, seq_len(qx$rank), drop = FALSE]; Pi <- tcrossprod(Qx)
  eps <- drop((diag(ic$M) - Pi) %*% ic$yk)
  Dc <- nb$D[ic$cand, ic$cand]; K <- pmax(1 - Dc / omega, 0)
  cc <- tcrossprod(ic$ck)
  hac <- sum(K * cc * tcrossprod(eps))
  R <- sum(K * cc * (diag(ic$M) - Pi))
  f <- hac / (ic$sigma2 * R)
  if (!is.finite(f) || R <= 0 || hac <= 0) NA_real_ else f
}

## ---- Sections 2-3: image data for every requested pair --------------------------------------------------------

## cells: imageID, cellType, x, y; pairs: data.frame(from, to). Returns, per pair, the image rows (O, E, n, v with
## v already multiplied by psi when label_clustering = TRUE), as the package's excess_image_data().
ref_image_data <- function(cells, pairs, r, label_clustering = TRUE, self_count = "spec", effect = "count") {
  imgs <- sort(unique(as.character(cells$imageID)))
  counts_fun <- if (effect == "allocation") ref_image_counts_alloc else ref_image_counts
  per_img <- lapply(imgs, function(im) {
    z <- cells[cells$imageID == im, ]; nb <- ref_neighbours(z$x, z$y, r)
    ic <- lapply(seq_len(nrow(pairs)), function(j) counts_fun(as.character(z$cellType), nb, pairs$from[j], pairs$to[j]))
    psi_t <- if (label_clustering) vapply(ic, ref_psi_tilde, 0, nb = nb, omega = 2 * r, self_count = self_count) else rep(NA_real_, nrow(pairs))
    # pooling: per target type, the median over the reference types of every requested pair with that target, floor 1
    psi <- vapply(seq_len(nrow(pairs)), function(j) { k <- pairs$to == pairs$to[j]; z <- psi_t[k]
      if (!label_clustering || all(is.na(z))) 1 else max(1, stats::median(z, na.rm = TRUE)) }, 0)
    list(ic = ic, psi = psi) })
  stats::setNames(lapply(seq_len(nrow(pairs)), function(j) {
    d <- do.call(rbind, lapply(seq_along(imgs), function(i) { ic <- per_img[[i]]$ic[[j]]
      data.frame(imageID = imgs[i], O = ic$O, E = ic$E, n = ic$n, v = max(ic$v, 0) * per_img[[i]]$psi[j]) }))
    d[d$n > 0 & d$v > 0, , drop = FALSE] }), paste0(pairs$from, "__", pairs$to))
}

## ---- Section 3: the frailty GEE, written as a general GEE ----------------------------------------------------

## Working covariance and precision of each unit: V_i = diag(v) + tau2 n n^T (Section 3.1).
ref_unit_V <- function(d, idx, tau2) lapply(idx, function(u) diag(d$v[u], length(u)) + tau2 * tcrossprod(d$n[u]))

## GEE fit for beta = (delta_1, delta_2) at a given tau2: beta = F^-1 sum D_i^T W_i Y_i.
ref_gee_fit <- function(d, unit, group, tau2) {
  idx <- split(seq_len(nrow(d)), factor(unit, levels = unique(unit)))
  g_u <- vapply(idx, function(u) group[u[1]], 0)
  V <- ref_unit_V(d, idx, tau2); W <- lapply(V, solve)
  Dm <- lapply(seq_along(idx), function(k) { u <- idx[[k]]; m <- matrix(0, length(u), 2); m[, g_u[k] + 1] <- d$n[u]; m })
  Y <- lapply(idx, function(u) d$O[u] - d$E[u])
  Fm <- Reduce(`+`, Map(function(D, W) t(D) %*% W %*% D, Dm, W))
  Bm <- solve(Fm)
  beta <- drop(Bm %*% Reduce(`+`, Map(function(D, W, y) t(D) %*% W %*% y, Dm, W, Y)))
  e <- Map(function(D, y) drop(y - D %*% beta), Dm, Y)
  list(idx = idx, g_u = g_u, V = V, W = W, D = Dm, Y = Y, F = Fm, B = Bm, beta = beta, e = e)
}

## Pearson statistic (eq. pearson, second form): sum_i w_i (delta~_i - delta^_g)^2, w_i = (1 / I_i + tau2)^-1, with
## delta~_i the information-weighted patient summary (eq. score) and delta^_g from the general GEE fit.
ref_pearson <- function(d, unit, group, tau2) {
  f <- ref_gee_fit(d, unit, group, tau2)
  sum(vapply(seq_along(f$idx), function(k) { u <- f$idx[[k]]
    I <- sum(d$n[u]^2 / d$v[u]); dt <- sum(d$n[u] * (d$O[u] - d$E[u]) / d$v[u]) / I
    (dt - f$beta[f$g_u[k] + 1])^2 / (1 / I + tau2) }, 0))
}

## Paule-Mandel (Section 3.2): the root of Q(tau2) = m - 2, or 0 if Q(0) <= m - 2.
ref_tau2 <- function(d, unit, group) {
  df <- length(unique(unit)) - 2; q0 <- ref_pearson(d, unit, group, 0)
  if (df < 1 || q0 <= df) return(0)
  up <- 1e-4; while (ref_pearson(d, unit, group, up) > df) up <- up * 4
  stats::uniroot(function(t) ref_pearson(d, unit, group, t) - df, c(0, up), tol = 1e-12)$root
}

msqrt <- function(S, pow) { e <- eigen((S + t(S)) / 2, symmetric = TRUE); e$vectors %*% (t(e$vectors) * e$values^pow) }

## CR2 (eq. A-whitened) and Satterthwaite df (Propositions 6-8) from their general definitions: explicit Phi_i,
## explicit stacked hat matrix H, loading vectors gamma_i, Omega = Gamma^T Sigma Gamma.
ref_cr2 <- function(f, ell = c(-1, 1)) {
  m <- length(f$idx); J <- sum(lengths(f$idx)); off <- cumsum(c(0, lengths(f$idx)))   # images stacked unit by unit
  Dall <- do.call(rbind, f$D); Wall <- ref_bdiag(f$W); Sig <- ref_bdiag(f$V)
  H <- Dall %*% f$B %*% t(Dall) %*% Wall
  Phi <- lapply(seq_len(m), function(i) { Wh <- msqrt(f$W[[i]], 0.5); Whi <- msqrt(f$W[[i]], -0.5)
    Whi %*% msqrt(diag(nrow(Wh)) - Wh %*% f$D[[i]] %*% f$B %*% t(f$D[[i]]) %*% Wh, -0.5) %*% Wh })
  zeta <- vapply(seq_len(m), function(i) sum(ell * (f$B %*% t(f$D[[i]]) %*% f$W[[i]] %*% Phi[[i]] %*% f$e[[i]])), 0)
  Gam <- vapply(seq_len(m), function(i) { r <- (off[i] + 1):off[i + 1]
    drop(t(ell) %*% f$B %*% t(f$D[[i]]) %*% f$W[[i]] %*% Phi[[i]] %*% (diag(J) - H)[r, , drop = FALSE]) }, numeric(J))
  Om <- t(Gam) %*% Sig %*% Gam
  list(V = sum(zeta^2), df = sum(diag(Om))^2 / sum(Om^2), Omega = Om)
}

## One pair: the test of delta_2 = delta_1 (Section 4.4), and the Hartung-Knapp option (Section 4.5).
ref_excess_test <- function(d, unit, group, variance = c("cr2", "hartung_knapp")) {
  variance <- match.arg(variance)
  tau2 <- ref_tau2(d, unit, group); f <- ref_gee_fit(d, unit, group, tau2); cr <- ref_cr2(f)
  est <- f$beta[2] - f$beta[1]; V <- cr$V; df <- cr$df
  if (variance == "hartung_knapp") {
    m <- length(f$idx); X <- ref_pearson(d, unit, group, tau2); Vm <- f$B[1, 1] + f$B[2, 2]
    V <- max(V, Vm * max(1, X / (m - 2))); df <- m - 2 }
  c(delta_ref = f$beta[1], delta_comp = f$beta[2], excess_difference = est, tau2 = tau2, se = sqrt(V), df = df,
    p_value = 2 * stats::pt(-abs(est / sqrt(V)), df))
}

## Whole dataset: cells with imageID, cellType, x, y, condition (two levels) and optionally subject.
ref_spicy_cell <- function(cells, pairs, r, label_clustering = TRUE, variance = "cr2", self_count = "spec",
                           effect = "count") {
  ims <- ref_image_data(cells, pairs, r, label_clustering, self_count, effect)
  info <- unique(cells[, intersect(c("imageID", "subject", "condition"), names(cells))])
  lev <- levels(factor(cells$condition))
  do.call(rbind, lapply(names(ims), function(k) { d <- ims[[k]]; i <- match(d$imageID, info$imageID)
    unit <- if (!is.null(info$subject)) as.character(info$subject[i]) else d$imageID
    group <- as.numeric(as.character(info$condition[i]) == lev[2])
    ft <- strsplit(k, "__", fixed = TRUE)[[1]]
    ok <- all(vapply(0:1, function(g) length(unique(unit[group == g])) >= 2, TRUE))
    res <- if (ok) ref_excess_test(d, unit, group, variance) else
      c(delta_ref = NA, delta_comp = NA, excess_difference = NA, tau2 = NA, se = NA, df = NA, p_value = NA)
    data.frame(from = ft[1], to = ft[2], t(res)) }))
}

ref_bdiag <- function(ms) {
  n <- vapply(ms, nrow, 1L); out <- matrix(0, sum(n), sum(n)); o <- cumsum(c(0, n))
  for (i in seq_along(ms)) out[(o[i] + 1):o[i + 1], (o[i] + 1):o[i + 1]] <- ms[[i]]
  out
}

## A general design (new_methods.pdf, Section 1) from the definitions: dense GLS with V_i = diag(v) + tau2 n n',
## the Bell-McCaffrey adjustment with symmetric square roots, and Omega = Gamma' Sigma Gamma.
ref_design <- function(d, unit, Z, cvec, tau2) {
  idx <- split(seq_len(nrow(d)), factor(unit, levels = unique(unit)))
  V <- lapply(idx, function(u) diag(d$v[u], length(u)) + tau2 * tcrossprod(d$n[u])); W <- lapply(V, solve)
  D <- lapply(idx, function(u) d$n[u] * Z[u, , drop = FALSE]); Y <- lapply(idx, function(u) d$O[u] - d$E[u])
  B <- solve(Reduce(`+`, Map(function(D, W) t(D) %*% W %*% D, D, W)))
  theta <- drop(B %*% Reduce(`+`, Map(function(D, W, y) t(D) %*% W %*% y, D, W, Y)))
  f <- list(idx = idx, V = V, W = W, D = D, B = B, e = Map(function(D, y) drop(y - D %*% theta), D, Y))
  cr <- ref_cr2(f, cvec)
  est <- sum(cvec * theta)
  c(estimate = est, se = sqrt(cr$V), df = cr$df, p = 2 * stats::pt(-abs(est / sqrt(cr$V)), cr$df))
}

## A small data set: n_pat patients per group, two images each, six types; in group B, a fraction of the
## T cells sit next to tumour cells (T cells recruited around tumour cells).
sim_cells <- function(seed = 1, n_pat = 6, n_cells = 350, side = 400) {
  set.seed(seed)
  types <- c("tumour", "T", "B", "macro", "stroma", "rare"); prob <- c(.3, .15, .15, .15, .2, .05)
  out <- list(); k <- 0
  for (g in c("A", "B")) for (p in seq_len(n_pat)) for (im in 1:2) {
    k <- k + 1; n <- n_cells + sample(-50:50, 1)
    z <- data.frame(x = runif(n, 0, side), y = runif(n, 0, side), cellType = sample(types, n, TRUE, prob))
    if (g == "B") { tum <- which(z$cellType == "tumour"); tc <- which(z$cellType == "T")
      mv <- tc[runif(length(tc)) < 0.4]; h <- sample(tum, length(mv), TRUE)
      z$x[mv] <- z$x[h] + runif(length(mv), -6, 6); z$y[mv] <- z$y[h] + runif(length(mv), -6, 6) }
    z$imageID <- sprintf("%s%02d_%d", g, p, im); z$patient <- sprintf("%s%02d", g, p); z$condition <- g
    out[[k]] <- z }
  d <- do.call(rbind, out); d$condition <- factor(d$condition); d
}
