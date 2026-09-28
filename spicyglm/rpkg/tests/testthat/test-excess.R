## effect = "excess" (math.tex, Section 18).

# images with a fraction `f` of B cells recruited to within 15 of A cells
make_excess_cells <- function(seed = 3, n_img = 12, f = c(0.1, 0.3), n = c(150, 150, 300, 400), side = NULL) {
  set.seed(seed)
  do.call(rbind, lapply(seq_len(n_img), function(i) {
    g <- 1 + (i > n_img / 2); w <- if (is.null(side)) runif(1, 500, 900) else side[i]
    tt <- rep(c("A", "B", "C", "D"), rpois(4, n))
    x <- runif(length(tt), 0, w); y <- runif(length(tt), 0, w)
    a <- which(tt == "A"); b <- which(tt == "B"); mv <- b[runif(length(b)) < f[g]]
    h <- a[sample.int(length(a), length(mv), TRUE)]; rr <- 15 * sqrt(runif(length(mv))); th <- runif(length(mv), 0, 2 * pi)
    x[mv] <- x[h] + rr * cos(th); y[mv] <- y[h] + rr * sin(th)
    data.frame(x = x, y = y, cellType = tt, imageID = sprintf("i%02d", i), patient = sprintf("p%02d", (i + 1) %/% 2),
               response = c("alpha", "beta")[g], stringsAsFactors = FALSE)
  }))
}

msym <- function(M, pw) { e <- eigen((M + t(M)) / 2, symmetric = TRUE); e$vectors %*% diag(e$values^pw, length(e$values)) %*% t(e$vectors) }

excess_ctx_image <- function(z, knn, r = 30, k = 10) {
  tl <- sort(unique(z$cellType))
  d <- dataset_create(z$x, z$y, match(z$cellType, tl) - 1L, c(0L, nrow(z)), length(tl))
  if (knn) dataset_build_knn(d, k, 1L) else dataset_build_radius_index(d, r)
  list(tl = tl, ctx = list(family = if (knn) "binomial" else "poisson",
                           totals = dataset_pair_neighbour_totals(d, knn, length(tl)),
                           out_sq_totals = dataset_pair_neighbour_out_sq_totals(d, knn, length(tl)),
                           counts = matrix(tabulate(match(z$cellType, tl), length(tl)), 1)))
}

# c_j for every cell: the number of `from` cells among j's neighbours
neighbour_counts <- function(z, from, knn, r = 30, k = 10) {
  dm <- as.matrix(stats::dist(z[, c("x", "y")])); diag(dm) <- Inf
  isA <- z$cellType == from
  if (knn) apply(dm, 1, function(dr) sum(isA[order(dr)[seq_len(k)]])) else rowSums(dm[, isA, drop = FALSE] <= r)
}

test_that("O, its random-labelling mean and variance match a direct computation and relabelling", {
  z <- make_excess_cells()[make_excess_cells()$imageID == "i02", ]
  for (knn in c(FALSE, TRUE)) {
    im <- excess_ctx_image(z, knn)
    for (pr in list(c("A", "B"), c("C", "A"), c("D", "C"))) {
      got <- spicyglm:::excess_image_data(im$ctx, match(pr[1], im$tl) - 1L, match(pr[2], im$tl) - 1L)
      cj <- neighbour_counts(z, pr[1], knn); cand <- z$cellType != pr[1]; y <- z$cellType == pr[2]
      cb <- cj[cand]; nB <- sum(y); M <- length(cb); p <- nB / M
      expect_equal(got$O, sum(cj[y]))
      expect_equal(got$E, p * sum(cb), tolerance = 1e-12)
      expect_equal(got$v, nB * (M - nB) / (M * (M - 1)) * sum((cb - mean(cb))^2), tolerance = 1e-10)
      set.seed(2); sims <- replicate(3000, sum(cb[sample(M, nB)]))
      expect_equal(mean(sims), got$E, tolerance = 0.02)
      expect_equal(var(sims), got$v, tolerance = 0.1)
    }
  }
})

test_that("recruiting a fraction f of TARGET cells gives an excess near f at any density", {
  # the same recruitment in sparse (large window) and dense (small window) images
  for (side in c(1500, 500)) {
    cells <- make_excess_cells(seed = 5, n_img = 12, f = c(0.2, 0.2), n = c(300, 300, 600, 800), side = rep(side, 12))
    out <- spicy_glm(cells, condition = "response", r = 30, effect = "excess", from = "A", to = "B", label_clustering = FALSE)
    expect_equal(out$results$coef_ref, 0.2, tolerance = 0.25)     # relative: 0.15 - 0.25
    expect_equal(out$results$coef_comp, 0.2, tolerance = 0.25)
  }
})

test_that("the closed-form CR2 equals the matrix CR2 of the linear GEE with several images per unit", {
  cells <- make_excess_cells(seed = 9, n_img = 16)
  out <- spicy_glm(cells, condition = "response", r = 30, subject = "patient", effect = "excess",
                   from = "A", to = "B", label_clustering = FALSE)
  # rebuild the image data and refit by matrices
  cs <- cells[order(cells$imageID), ]; labs <- sort(unique(cs$imageID)); tl <- unique(cells$cellType)
  img <- match(cs$imageID, labs) - 1L; off <- as.integer(c(0, cumsum(tabulate(img + 1L, nbins = length(labs)))))
  d <- dataset_create(cs$x, cs$y, match(cs$cellType, tl) - 1L, off, length(tl)); dataset_build_radius_index(d, 30)
  ctx <- list(family = "poisson", totals = dataset_pair_neighbour_totals(d, FALSE, length(tl)),
              out_sq_totals = dataset_pair_neighbour_out_sq_totals(d, FALSE, length(tl)),
              counts = unclass(table(factor(img, levels = seq_along(labs) - 1L), factor(match(cs$cellType, tl) - 1L, levels = seq_along(tl) - 1L))))
  im <- spicyglm:::excess_image_data(ctx, match("A", tl) - 1L, match("B", tl) - 1L)
  unit <- tapply(cs$patient, img, `[`, 1)[im$img + 1L]; grp <- tapply(cs$response, img, `[`, 1)[im$img + 1L] == "beta"
  tau2 <- out$frailty$tau2
  V <- 0; EV <- 0; trsq <- 0; est <- numeric(2)
  for (g in c(FALSE, TRUE)) {
    k <- grp == g; x <- im$n[k]; y <- im$O[k] - im$E[k]; u <- unit[k]; us <- unique(u)
    W <- lapply(us, function(s) solve(diag(im$v[k][u == s], sum(u == s)) + tau2 * tcrossprod(x[u == s])))
    B <- sum(vapply(seq_along(us), function(j) { xs <- x[u == us[j]]; drop(t(xs) %*% W[[j]] %*% xs) }, 0))
    delta <- sum(vapply(seq_along(us), function(j) { s <- u == us[j]; drop(t(x[s]) %*% W[[j]] %*% y[s]) }, 0)) / B
    est[g + 1] <- delta
    for (j in seq_along(us)) { s <- u == us[j]; xs <- x[s]; e <- y[s] - delta * xs
      # Bell-McCaffrey: A = W^-1/2 (I - W^1/2 x B^-1 x' W^1/2)^-1/2 W^1/2
      Wh <- msym(W[[j]], 0.5); Wmh <- msym(W[[j]], -0.5)
      A <- Wmh %*% msym(diag(sum(s)) - Wh %*% tcrossprod(xs) %*% Wh / B, -0.5) %*% Wh
      V <- V + drop(t(xs) %*% W[[j]] %*% A %*% e)^2 / B^2 }
  }
  expect_equal(unname(out$results$coef_ref), est[1], tolerance = 1e-8)
  expect_equal(unname(out$results$coef_comp), est[2], tolerance = 1e-8)
  expect_equal(out$frailty$se^2, V, tolerance = 1e-6)
})

test_that("kNN excess, subjects and argument checks", {
  cells <- make_excess_cells()
  out <- spicy_glm(cells, condition = "response", family = "binomial", k = 10, effect = "excess")
  expect_true(all(c("excess_difference", "coef_ref", "coef_comp", "p_value") %in% names(out$results)))
  expect_equal(nrow(out$results), 16L)                      # every ordered pair, self-pairs included
  expect_true(all(is.finite(out$results$p_value)))
  ab <- out$results[out$results$from == "A" & out$results$to == "B", ]
  expect_gt(ab$excess_difference, 0)
  expect_error(spicy_glm(cells, condition = "response", r = 30, effect = "excess", sigma = 100), "random-labelling")
  expect_error(spicy_glm(cells, condition = "response", r = 30, effect = "excess", test = "moderated"), "cell-level")
  indep <- spicy_glm(cells, condition = "response", r = 30, effect = "excess", frailty = FALSE, from = "A", to = "B")
  expect_equal(indep$frailty$tau2, 0)
})

test_that("moderation keeps each pair's estimate and adds prior df", {
  cells <- make_excess_cells(seed = 4, n_img = 12, n = c(150, 150, 300, 400))
  a <- spicy_glm(cells, condition = "response", r = 30, effect = "excess")
  b <- spicy_glm(cells, condition = "response", r = 30, effect = "excess", moderate = TRUE)
  m <- match(paste(a$results$from, a$results$to), paste(b$results$from, b$results$to))
  expect_equal(b$results$excess_difference[m], a$results$excess_difference, tolerance = 1e-12)
  fa <- a$frailty[order(a$frailty$from, a$frailty$to), ]; fb <- b$frailty[order(b$frailty$from, b$frailty$to), ]
  expect_true(all(fb$df >= fa$df - 1e-8))
  pr <- attr(b$frailty, "prior"); expect_true(all(c("tau2", "d0", "s0") %in% names(pr)))
  expect_true(all(is.finite(b$results$p_value)))
})
