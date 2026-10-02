## null = "random_labelling" and frailty = TRUE (math.tex, Section 17).

make_frailty_cells <- function(seed = 8, n_img = 12) {
  set.seed(seed)
  do.call(rbind, lapply(seq_len(n_img), function(i) {
    n <- rpois(4, c(220, 160, 90, 50)); tt <- rep(c("A", "B", "C", "D"), n)
    width <- runif(1, 400, 1000)                       # tissue fills a varying share of the frame
    left <- runif(length(tt)) < ifelse(tt %in% c("A", "B"), 0.35 + 0.1 * (i > n_img / 2), 0.6)
    data.frame(x = ifelse(left, runif(length(tt), 0, width / 2), runif(length(tt), width / 2, width)),
               y = runif(length(tt), 0, 600), cellType = tt, imageID = sprintf("i%02d", i),
               patient = sprintf("p%02d", (i + 1) %/% 2),
               response = if (i <= n_img / 2) "alpha" else "beta", stringsAsFactors = FALSE)
  }))
}

one_image <- function(cells, id, r) {
  z <- cells[cells$imageID == id, ]; tl <- sort(unique(cells$cellType))
  d <- dataset_create(z$x, z$y, match(z$cellType, tl) - 1L, c(0L, nrow(z)), length(tl))
  dataset_build_radius_index(d, r)
  list(z = z, tl = tl, d = d, dm = { m <- as.matrix(stats::dist(z[, c("x", "y")])); diag(m) <- Inf; m })
}

test_that("the random-labelling offset is each REF cell's candidate neighbours times the TARGET share", {
  cells <- make_frailty_cells(); im <- one_image(cells, "i01", 40)
  for (pr in list(c("A", "B"), c("C", "A"), c("B", "B"))) {
    f <- match(pr[1], im$tl) - 1L; t <- match(pr[2], im$tl) - 1L
    md <- dataset_rl_model_data(im$d, f, t)
    a <- which(im$z$cellType == pr[1]); close <- im$dm[a, , drop = FALSE] <= 40
    self <- pr[1] == pr[2]
    cand <- if (self) rep(TRUE, nrow(im$z)) else im$z$cellType != pr[1]
    share <- if (self) (length(a) - 1) / (nrow(im$z) - 1) else sum(im$z$cellType == pr[2]) / sum(cand)
    expect_equal(md$n, unname(as.integer(rowSums(close[, im$z$cellType == pr[2], drop = FALSE]))))
    expect_equal(md$density, unname(share * rowSums(close[, cand, drop = FALSE])), tolerance = 1e-12)
  }
})

test_that("phi is the random-labelling variance of the count over its mean", {
  cells <- make_frailty_cells(); im <- one_image(cells, "i03", 40)
  T <- length(im$tl)
  tot <- dataset_pair_neighbour_totals(im$d, FALSE, T); sq <- dataset_pair_neighbour_sq_totals(im$d, FALSE, T)
  counts <- matrix(tabulate(match(im$z$cellType, im$tl), T), 1)
  f <- match("A", im$tl); t <- match("C", im$tl)
  phi <- spicyR:::frailty_phi(tot, sq, counts, f, t, FALSE)
  # exact moments over all relabellings of the non-A cells with n_C fixed
  isA <- im$z$cellType == "A"; cb <- colSums(im$dm[isA, , drop = FALSE] <= 40)[!isA]
  nC <- sum(im$z$cellType == "C"); N <- length(cb); p <- nC / N
  v <- nC * (N - nC) / (N * (N - 1)) * sum((cb - mean(cb))^2)
  expect_equal(phi, v / (p * sum(cb)), tolerance = 1e-12)
  set.seed(1); sims <- replicate(4000, sum(cb[sample(N, nC)]))
  expect_equal(var(sims) / mean(sims), phi, tolerance = 0.1)
})

test_that("with tau2 = 0 and phi = 1 the frailty CR2 is the cell-level closed form", {
  cells <- make_frailty_cells()
  cs <- cells[order(cells$imageID), ]; labs <- sort(unique(cs$imageID)); img <- match(cs$imageID, labs) - 1L
  tl <- sort(unique(cs$cellType)); off <- as.integer(c(0, cumsum(tabulate(img + 1L, nbins = length(labs)))))
  d <- dataset_create(cs$x, cs$y, match(cs$cellType, tl) - 1L, off, length(tl)); dataset_build_radius_index(d, 40)
  area <- dataset_image_areas(d, "convex"); grp <- as.integer(tapply(cs$response, img, `[`, 1) == "beta")
  md <- dataset_poisson_model_data(d, area, 0L, 1L)
  ref <- fit_pair_poisson_cpp(md$image, md$image, grp[md$image + 1L], md$n, md$density, "firth", "fast", FALSE)
  im <- data.frame(img = sort(unique(md$image)), O = as.numeric(tapply(md$n, md$image, sum)),
                   E = as.numeric(tapply(md$density, md$image, sum)), phi = 1)
  ctx <- list(family = "poisson", image_cluster = seq_along(labs) - 1L, image_group = grp)
  unit <- factor(ctx$image_cluster[im$img + 1L]); g <- grp[im$img + 1L]
  fits <- lapply(0:1, function(k) spicyR:::frailty_fit_group(im[g == k, ], droplevels(unit[g == k]), 0, FALSE))
  cr <- lapply(fits, spicyR:::frailty_cr2)
  expect_equal(c(fits[[1]]$beta, fits[[2]]$beta), ref$beta, tolerance = 1e-10)
  expect_equal(cr[[1]]$V + cr[[2]]$V, ref$v_hat, tolerance = 1e-10)
  expect_equal((cr[[1]]$EV + cr[[2]]$EV)^2 / (cr[[1]]$trsq + cr[[2]]$trsq), ref$df, tolerance = 1e-10)
})

test_that("frailty = TRUE runs for the Poisson (both nulls) and binomial designs, with subjects", {
  cells <- make_frailty_cells()
  for (args in list(list(r = 40), list(r = 40, null = "random_labelling"), list(family = "binomial", k = 10))) {
    out <- do.call(spicy_glm, c(list(cells, condition = "response", subject = "patient", frailty = TRUE, effect = "ratio"), args))
    expect_true(nrow(out$results) > 0)
    expect_true(all(is.finite(out$results$p_value)))
    expect_true(all(out$frailty$tau2 >= 0) && all(out$frailty$df > 0))
    expect_setequal(paste(out$frailty$from, out$frailty$to), paste(out$results$from, out$results$to))
  }
  rl <- spicy_glm(effect = "ratio", cells, condition = "response", r = 40, null = "random_labelling")
  expect_equal(nrow(rl$results) + nrow(rl$skipped), 16)           # directional pairs
  expect_error(spicy_glm(effect = "ratio", cells, condition = "response", family = "binomial", k = 10, null = "random_labelling"),
               "applies to family")
})

test_that("moderate = TRUE shrinks the CR2 variances toward the frailty-model variance", {
  cells <- make_frailty_cells()
  plain <- spicy_glm(effect = "ratio", cells, condition = "response", r = 40, null = "random_labelling", frailty = TRUE)
  mod <- spicy_glm(effect = "ratio", cells, condition = "response", r = 40, null = "random_labelling", frailty = TRUE, moderate = TRUE)
  pr <- attr(mod$frailty, "prior")
  expect_named(pr, c("tau2", "d0", "s0"))
  expect_equal(unname(pr["tau2"]), stats::median(plain$frailty$tau2))
  key <- paste(mod$frailty$from, mod$frailty$to); kp <- match(key, paste(plain$frailty$from, plain$frailty$to))
  if (is.finite(pr["d0"])) expect_true(all(mod$frailty$df >= plain$frailty$df[kp] - 1e-8))
  expect_error(spicy_glm(effect = "ratio", cells, condition = "response", r = 40, moderate = TRUE), "needs frailty")
})

test_that("the Kontextual phi is the variance over relabellings of the context", {
  set.seed(5); n <- 600
  z <- data.frame(x = runif(n, 0, 500), y = runif(n, 0, 500))
  z$t <- ifelse(z$x < 200, sample(c("A", "B", "C", "D"), n, TRUE, c(.3, .3, .2, .2)),
                sample(c("A", "B", "C", "D"), n, TRUE, c(.15, .15, .35, .35)))
  tl <- c("A", "B", "C", "D"); ctxt <- match(c("B", "C", "D"), tl) - 1L
  total <- function(types) { d <- dataset_create(z$x, z$y, match(types, tl) - 1L, c(0L, n), 4L)
    dataset_build_radius_index(d, 40); dataset_build_context(d, ctxt, "convex", FALSE)
    list(d = d, O = sum(dataset_kontextual_model_data(d, 0L, 1L)$n)) }
  r0 <- total(z$t); S <- dataset_weighted_phi_sums(r0$d, 0L, 1L, 0L, FALSE)
  p <- sum(z$t == "B") / S[3, 1]; v <- p * (1 - p) * S[3, 1] / (S[3, 1] - 1) * (S[2, 1] - S[1, 1]^2 / S[3, 1])
  cand <- which(z$t %in% c("B", "C", "D")); lab <- z$t[cand]
  set.seed(2); sims <- replicate(800, { tt <- z$t; tt[cand] <- sample(lab); total(tt)$O })
  expect_equal(mean(sims), p * S[1, 1], tolerance = 0.02)
  expect_equal(var(sims), v, tolerance = 0.15)
})

test_that("frailty = TRUE runs for the Kontextual and inhomogeneous designs", {
  cells <- make_frailty_cells()
  for (args in list(list(r = 40, parent = c("B", "C", "D")), list(family = "binomial", k = 10, parent = c("B", "C", "D")),
                    list(r = 40, sigma = 120))) {
    out <- do.call(spicy_glm, c(list(cells, condition = "response", subject = "patient", frailty = TRUE), args))
    expect_true(nrow(out$results) > 0)
    expect_true(all(is.finite(out$results$p_value)))
  }
})
