## Adjusted designs with a column that varies between a patient's images (extension A/B; DECISIONS 4 Nov 2026 17:10
## and 6 Nov 14:37): the image variance by Paule-Mandel within patients after the covariates (B2), the patient mean of
## each such column (Mundlak), and tau2 by Paule-Mandel on the patient summaries. The core against the plain-R
## definitions (ref_adjusted_design in helper-reference.R) and the prototype, and the unadjusted test unchanged.

## Row-level data: 24 patients with 1 to 6 images, a covariate x that varies within patients and is shifted by
## condition, s2 the rows' sampling variances (the prototype's simulation design, analysis/explore/2026-11_mundlak).
ad_data <- function() {
  set.seed(20261106); m <- 24; J <- rep(c(1, 3, 6, 2, 4, 5), 4); id <- rep(seq_len(m), J); N <- sum(J)
  cond <- rep(0:1, each = m / 2); u <- rnorm(m, 0.5 * cond)
  x <- sqrt(0.5) * u[id] + sqrt(0.5) * rnorm(N); s2 <- exp(rnorm(N, log(0.01), 0.9))
  y <- 0.05 * cond[id] + 0.05 * x + 0.1 * u[id] + rnorm(m, 0, 0.1)[id] + rnorm(N, 0, sqrt(0.005)) + rnorm(N, 0, sqrt(s2))
  list(raw = data.frame(img = 0:(N - 1), O = y, E = 0, n = 1, v = s2), id = id, x = x, cond = cond, m = m, N = N)
}

test_that("adjusted design: the core gives the prototype's numbers", {
  # analysis/explore/2026-11_mundlak/mundlak_lib.R (omega2_B2, mundlak_cols(weighted = FALSE), tau2_between_rank) and
  # 2026-11_tau2_adjusted/tau2_lib.R (cr2_test) on this data set, 11 Nov 2026
  d <- ad_data()
  h <- spicyR:::stats_image_heterogeneity_rows(d$raw, d$id - 1L)
  Z <- cbind(d$cond[d$id] == 0, d$cond[d$id] == 1, d$x - mean(d$x)) * 1
  a <- spicyR:::stats_adjusted_design(h, d$id - 1L, Z, 2L, TRUE)
  expect_true(a$image_level)
  expect_identical(a$mean_of, 3L)
  r <- spicyR:::stats_design_tests(a$rows, d$id - 1L, d$m, a$Z, rbind(c(-1, 1, 0, 0), c(0, 0, 1, 0)), a$tau2, FALSE, 0L)
  expect_equal(c(a$omega2, a$tau2), c(0.00359010957756, 0.00897012977309), tolerance = 1e-8)
  expect_equal(c(r[[1]]$estimate, r[[1]]$se, r[[1]]$df, r[[1]]$p), c(0.131356981759, 0.046972032999, 17.7312613517,
                                                                     0.0120465213583), tolerance = 1e-8)
  expect_equal(r[[2]]$estimate, 0.00285477533688, tolerance = 1e-8)   # the within-patient slope of x
})

test_that("adjusted design: omega2, the mean columns and tau2 match the definitions", {
  d <- ad_data()
  h <- spicyR:::stats_image_heterogeneity_rows(d$raw, d$id - 1L)
  expect_equal(h$s2, d$raw$v / d$raw$n^2)
  set.seed(2); x2 <- rnorm(d$m)[d$id]; x3 <- rnorm(d$N)   # a patient-level column (no mean) and a second image-level one
  Z <- cbind(d$cond[d$id] == 0, d$cond[d$id] == 1, d$x - mean(d$x), x2 - mean(x2), x3 - mean(x3)) * 1
  a <- spicyR:::stats_adjusted_design(h, d$id - 1L, Z, 2L, TRUE)
  ref <- ref_adjusted_design(d$raw, d$id, Z, 2)
  expect_identical(a$mean_of, c(3L, 5L))
  expect_equal(unname(a$Z), unname(ref$Z), tolerance = 1e-12)
  expect_equal(c(a$omega2, a$tau2), c(ref$omega2, ref$tau2), tolerance = 1e-9)
  expect_equal(a$rows$v, ref$rows$v, tolerance = 1e-12)
  expect_equal(a$rows$O, ref$rows$O, tolerance = 1e-12)
  cv <- c(-1, 1, 0, 0, 0, 0, 0)
  t <- spicyR:::stats_design_tests(a$rows, d$id - 1L, d$m, a$Z, rbind(cv), a$tau2, FALSE, 0L)[[1]]
  expect_equal(c(t$estimate, t$se, t$df, t$p), unname(ref_design(ref$rows, d$id, ref$Z, cv, ref$tau2)), tolerance = 1e-8)
  # B2 differs from the unadjusted omega2 (x is absorbed there), and the unadjusted one is unchanged
  expect_gt(abs(a$omega2 - attr(h, "tau2_image")), 1e-6)
  expect_equal(attr(h, "tau2_image"), attr(ref_image_heterogeneity(d$raw, d$id), "tau2_image"), tolerance = 1e-12)
  # without the frailty tau2 is 0; rows without the image random effect (no s2) have no image variance
  expect_identical(spicyR:::stats_adjusted_design(h, d$id - 1L, Z, 2L, FALSE)$tau2, 0)
  b <- spicyR:::stats_adjusted_design(d$raw, d$id - 1L, Z, 2L, TRUE)
  expect_identical(b$omega2, 0)
  expect_equal(b$rows$v, d$raw$v)
  # every added column constant within patients: nothing is changed (the patient-level design, as before)
  expect_false(spicyR:::stats_adjusted_design(h, d$id - 1L, Z[, c(1, 2, 4)], 2L, TRUE)$image_level)
})

test_that("adjusted design: collinear image-level columns count once in both targets", {
  d <- ad_data()
  h <- spicyR:::stats_image_heterogeneity_rows(d$raw, d$id - 1L)
  # x and 2x span one direction: q_w = 1 and rank(zbar) is one less than the columns of the design
  Z <- cbind(d$cond[d$id] == 0, d$cond[d$id] == 1, d$x, 2 * d$x) * 1
  a <- spicyR:::stats_adjusted_design(h, d$id - 1L, Z, 2L, TRUE)
  ref <- ref_adjusted_design(d$raw, d$id, Z, 2)
  expect_equal(c(a$omega2, a$tau2), c(ref$omega2, ref$tau2), tolerance = 1e-9)
  one <- spicyR:::stats_adjusted_design(h, d$id - 1L, Z[, 1:3], 2L, TRUE)
  expect_equal(c(a$omega2, a$tau2), c(one$omega2, one$tau2), tolerance = 1e-9)
})

## spicy() on cells: two images per patient, so the abundance share and an image covariate vary within patients
cells_ad <- sim_cells()
set.seed(7); imz <- setNames(rnorm(length(unique(cells_ad$imageID))), unique(cells_ad$imageID))
cells_ad$z <- imz[cells_ad$imageID]
pat_ad <- unique(cells_ad$patient)
cells_ad$age <- setNames(seq(31, by = 3.5, length.out = length(pat_ad)), pat_ad)[cells_ad$patient]
cells_ad$cond3 <- factor(c("a", "b", "c")[(match(cells_ad$patient, pat_ad) %% 3) + 1])

## The pair's rows before the image random effect, the condition block (G columns) and the added columns, centred as
## the package centres them (covariates over the images, then every added column over the pair's rows).
ad_pair <- function(cond, f, t, cov = "z") {
  ctx <- spicyR:::.cell_context(cells_ad, cond, "patient", "imageID", "cellType", c("x", "y"))
  g <- spicyR:::.cell_graph(ctx, list(c(f, t)), r = 30, label_clustering = FALSE)
  raw <- spicyR:::.cell_rows(ctx, g, f, t)
  B <- spicyR:::.cell_condition_block(ctx, raw)
  info <- unique(cells_ad[, c("imageID", cov)]); zi <- info[[cov]][match(ctx$image_labels, info$imageID)]
  X <- cbind(abundance = spicyR:::.cell_share(ctx, raw, f), z = (zi - mean(zi))[raw$img + 1L])
  X <- sweep(X, 2, colMeans(X))
  list(raw = raw, Z = cbind(B, X), G = ncol(B), ctx = ctx)
}

test_that("spicy(), two conditions: the adjusted test (abundance and an image covariate) matches the definitions", {
  p <- ad_pair("condition", "tumour", "T")
  ref <- ref_adjusted_design(p$raw, p$raw$unit, p$Z, 2)
  cv <- c(-1, 1, rep(0, ncol(ref$Z) - 2))
  b <- ref_design(ref$rows, p$raw$unit, ref$Z, cv, ref$tau2)
  w <- ref_design(ref$rows, p$raw$unit, ref$Z, replace(0 * cv, 3, 1), ref$tau2)   # the abundance's within slope
  res <- suppressMessages(spicy(cells_ad, "condition", subject = "patient", r = 30, from = "tumour", to = "T",
                                covariates = "z", adjustAbundance = TRUE))$cellResults
  expect_identical(res$adjusted_for, "abundance+covariates")
  expect_equal(c(res$excess_difference, res$se, res$df, res$p_value), unname(b), tolerance = 1e-8)
  expect_equal(c(res$tau2, res$tau2_image), c(ref$tau2, ref$omega2), tolerance = 1e-9)
  expect_equal(c(res$abundance_effect, res$abundance_p_value), unname(w[c(1, 4)]), tolerance = 1e-8)
  # the core's every-pair path and the R path of one pair give the same test
  g <- spicyR:::.cell_graph(p$ctx, list(c("tumour", "T")), r = 30, label_clustering = FALSE); g$image_heterogeneity <- TRUE
  zx <- matrix(imz[p$ctx$image_labels] - mean(imz), ncol = 1, dimnames = list(NULL, "z"))
  o <- spicyR:::.cell_pair_test(p$ctx, g, "tumour", "T", TRUE, "cr2", TRUE, zx)
  expect_equal(c(o$test$difference, o$test$se, o$test$df, o$test$p, o$test$tau2, o$tau2_image),
               c(res$excess_difference, res$se, res$df, res$p_value, res$tau2, res$tau2_image), tolerance = 1e-12)
  # several radii: each radius's adjusted test is that radius's single-radius test
  rr <- suppressMessages(spicy(cells_ad, "condition", subject = "patient", r = c(20, 30), from = "tumour", to = "T",
                               covariates = "z", adjustAbundance = TRUE))$radiusResults
  expect_equal(rr$p_value[rr$r == 30], res$p_value, tolerance = 1e-12)
  expect_equal(rr$tau2_image[rr$r == 30], res$tau2_image, tolerance = 1e-12)
})

test_that("spicy(), three conditions and a continuous condition: the adjusted tests match the definitions", {
  p <- ad_pair("cond3", "tumour", "T")
  ref <- ref_adjusted_design(p$raw, p$raw$unit, p$Z, 3)
  q <- ncol(ref$Z)
  b <- ref_design(ref$rows, p$raw$unit, ref$Z, c(-1, 1, rep(0, q - 2)), ref$tau2)
  C <- rbind(c(-1, 1, rep(0, q - 2)), c(-1, 0, 1, rep(0, q - 3)))
  j <- ref_joint(ref$rows, p$raw$unit, ref$Z, C, ref$tau2)
  res <- suppressMessages(spicy(cells_ad, "cond3", subject = "patient", r = 30, from = "tumour", to = "T",
                                covariates = "z", adjustAbundance = TRUE))$cellResults
  x <- res[res$level == "b", ]
  expect_equal(c(x$excess_difference, x$se, x$df, x$p_value), unname(b), tolerance = 1e-8)
  expect_equal(c(x$overall_F, x$overall_df2, x$overall_p_value), unname(j), tolerance = 1e-8)
  expect_equal(c(x$tau2, x$tau2_image), c(ref$tau2, ref$omega2), tolerance = 1e-9)

  p <- ad_pair("age", "tumour", "T")
  ref <- ref_adjusted_design(p$raw, p$raw$unit, p$Z, 2)
  b <- ref_design(ref$rows, p$raw$unit, ref$Z, c(0, 1, rep(0, ncol(ref$Z) - 2)), ref$tau2)
  res <- suppressMessages(spicy(cells_ad, "age", subject = "patient", r = 30, from = "tumour", to = "T",
                                covariates = "z", adjustAbundance = TRUE))$cellResults
  expect_equal(c(res$excess_slope, res$se, res$df, res$p_value), unname(b), tolerance = 1e-8)
  expect_equal(c(res$tau2, res$tau2_image), c(ref$tau2, ref$omega2), tolerance = 1e-9)
})

test_that("without covariates, and with patient-level covariates only, results are those of spicyR 1.99.56", {
  # reference_1.99.56_adjusted.rds: spicy() of spicyR 1.99.56 on sim_cells() (patient-level age); bit-identical on x86-64 Linux,
  # compared to 1e-12 for other platforms' rounding
  ref <- readRDS(test_path("reference_1.99.56_adjusted.rds"))
  cells <- sim_cells()
  pat <- unique(cells$patient); cells$age <- setNames(seq(31, by = 3.5, length.out = length(pat)), pat)[cells$patient]
  s <- function(...) suppressMessages(spicy(cells, subject = "patient", ...))$cellResults
  # several radii: combine = "maxT", the default up to 1.99.57, as in the reference
  now <- list(allocation = s("condition", r = 30), count = s("condition", r = 30, effect = "count"),
              radii = s("condition", r = c(15, 30), combine = "maxT"), age = s("condition", r = 30, covariates = "age"),
              continuous = s("age", r = 30))
  for (k in names(ref)) {
    x <- now[[k]][, names(ref[[k]])]
    expect_equal(x, ref[[k]], tolerance = 1e-12, label = k)
  }
})
