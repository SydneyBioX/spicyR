## test = "moderated": image-level log O/E (or log odds ratio) summaries, compared
## between conditions with limma's trended empirical-Bayes moderated t.

make_mod_cells <- function(seed = 3, n_img = 12) {
  set.seed(seed)
  do.call(rbind, lapply(seq_len(n_img), function(i) {
    n <- rpois(6, c(200, 150, 120, 90, 60, 40)); tt <- rep(c("A", "B", "C", "D", "E", "F"), n)
    left <- runif(length(tt)) < ifelse(tt %in% c("A", "B"), 0.3 + 0.1 * (i > n_img / 2), 0.6)
    data.frame(x = ifelse(left, runif(length(tt), 0, 500), runif(length(tt), 500, 1000)),
               y = runif(length(tt), 0, 1000), cellType = tt, imageID = sprintf("i%02d", i),
               patient = sprintf("p%02d", (i + 1) %/% 2),
               response = if (i <= n_img / 2) "alpha" else "beta", stringsAsFactors = FALSE)
  }))
}

# the image summaries the moderated test works from, recomputed from the model data
reference_theta <- function(cells, family, r = NULL, k = NULL) {
  out <- spicy_glm(effect = "ratio", cells, condition = "response", r = r, k = k, family = family)
  pairs <- Map(c, out$results$from, out$results$to)
  cs <- cells[order(cells$imageID), ]
  labs <- sort(unique(cs$imageID)); img <- match(cs$imageID, labs) - 1L
  tl <- unique(cells$cellType); tc <- match(cs$cellType, tl) - 1L
  off <- as.integer(c(0, cumsum(tabulate(img + 1L, nbins = length(labs)))))
  d <- dataset_create(cs$x, cs$y, tc, off, length(tl))
  area <- dataset_image_areas(d, "convex")
  if (family == "poisson") dataset_build_radius_index(d, r) else dataset_build_knn(d, k)
  th <- t(vapply(pairs, function(p) {
    f <- match(p[1], tl) - 1L; t <- match(p[2], tl) - 1L
    if (family == "poisson") {
      md <- dataset_poisson_model_data(d, area, f, t)
      O <- tapply(md$n, factor(md$image, 0:(length(labs) - 1)), sum)
      log((O + 0.5) / tapply(md$density, factor(md$image, 0:(length(labs) - 1)), sum))
    } else {
      md <- dataset_binomial_model_data(d, f, t)
      im <- factor(md$image, 0:(length(labs) - 1))
      O <- tapply(md$n, im, sum); Tr <- tapply(rep(k, length(md$n)), im, sum)
      stats::qlogis((O + 0.5) / (Tr + 1)) - stats::qlogis(tapply(k * md$p0, im, sum) / Tr)
    }
  }, numeric(length(labs))))
  O <- t(vapply(pairs, function(p) {
    md <- if (family == "poisson") dataset_poisson_model_data(d, area, match(p[1], tl) - 1L, match(p[2], tl) - 1L)
          else dataset_binomial_model_data(d, match(p[1], tl) - 1L, match(p[2], tl) - 1L)
    as.numeric(tapply(md$n, factor(md$image, 0:(length(labs) - 1)), sum))
  }, numeric(length(labs))))
  list(pairs = pairs, theta = th, O = O,
       group = as.integer(tapply(cs$response, img, `[`, 1) == "beta"))
}

for (family in c("poisson", "binomial")) {
  test_that(paste("the", family, "moderated test is limma's eBayes(trend = TRUE)"), {
    skip_if_not_installed("limma")
    cells <- make_mod_cells()
    r <- if (family == "poisson") 50; k <- if (family == "binomial") 10
    ref <- reference_theta(cells, family, r, k)
    mod <- spicy_glm(effect = "ratio", cells, condition = "response", r = r, k = k, family = family, test = "moderated")
    A <- rowMeans(log1p(ref$O))
    e <- limma::eBayes(limma::lmFit(ref$theta, stats::model.matrix(~ ref$group)), trend = A, legacy = TRUE)
    key <- paste(mod$results$from, mod$results$to)
    j <- match(key, vapply(ref$pairs, paste, "", collapse = " "))
    expect_equal(mod$results$p_value, unname(e$p.value[j, 2]), tolerance = 1e-10)
    expect_equal(mod$results[[if (family == "poisson") "log_rate_ratio" else "log_odds_ratio"]],
                 unname(e$coefficients[j, 2]), tolerance = 1e-12)
    expect_equal(mod$moderation$prior_df, e$df.prior[1], tolerance = 1e-10)
    expect_true(all(mod$results$estimator == "moderated"))
  })
}

test_that("subjects are the units: image summaries are averaged within a subject", {
  skip_if_not_installed("limma")
  cells <- make_mod_cells()
  ref <- reference_theta(cells, "poisson", r = 50)
  mod <- spicy_glm(effect = "ratio", cells, condition = "response", r = 50, subject = "patient", test = "moderated")
  unit <- rep(1:6, each = 2)
  U <- t(apply(ref$theta, 1, function(z) tapply(z, unit, mean)))
  e <- limma::eBayes(limma::lmFit(U, stats::model.matrix(~ c(0, 0, 0, 1, 1, 1))),
                     trend = rowMeans(log1p(ref$O)), legacy = TRUE)
  j <- match(paste(mod$results$from, mod$results$to), vapply(ref$pairs, paste, "", collapse = " "))
  expect_equal(mod$results$p_value, unname(e$p.value[j, 2]), tolerance = 1e-10)
})

test_that("density adjustment removes a shared abundance slope and is invariant to condition shifts", {
  cells <- make_mod_cells()
  plain <- spicy_glm(effect = "ratio", cells, condition = "response", r = 50, test = "moderated")
  adj <- spicy_glm(effect = "ratio", cells, condition = "response", r = 50, test = "moderated", density_adjust = TRUE)
  expect_length(adj$moderation$density_slopes, 2)
  expect_null(plain$moderation$density_slopes)
  expect_setequal(paste(adj$results$from, adj$results$to), paste(plain$results$from, plain$results$to))

  # slopes are estimated within condition: adding a constant to every summary
  # of one condition leaves them unchanged
  s <- list(list(O = rep(10, 6), theta = c(1, 2, 3, 1.5, 2.5, 3.5)),
            list(O = rep(10, 6), theta = c(0.5, 1, 1.5, 0.2, 0.9, 1.1)))
  counts <- cbind(c(10, 20, 40, 15, 30, 60), c(5, 5, 5, 8, 8, 8))
  g <- c(0, 0, 0, 1, 1, 1)
  b <- spicyglm:::density_slopes(s, list(c(1, 2), c(2, 1)), counts, g)
  s2 <- lapply(s, function(z) { z$theta <- z$theta + 5 * g; z })
  expect_equal(spicyglm:::density_slopes(s2, list(c(1, 2), c(2, 1)), counts, g), b)
})

test_that("pairs with fewer than two units in a condition are skipped", {
  cells <- make_mod_cells()
  cells <- cells[!(cells$cellType == "F" & cells$response == "beta" & cells$imageID != "i12"), ]
  mod <- spicy_glm(effect = "ratio", cells, condition = "response", r = 50, test = "moderated")
  expect_true(all(mod$skipped$reason[grepl("F", paste(mod$skipped$from, mod$skipped$to))] == "one_patient_per_group"))
  expect_false(any(mod$results$from == "F" | mod$results$to == "F"))
})

test_that("Kontextual and inhomogeneous designs run through the per-pair route", {
  cells <- make_mod_cells()
  k <- spicy_glm(effect = "ratio", cells, condition = "response", family = "binomial", k = 10, parent = c("C", "D", "E"),
                 test = "moderated")
  expect_true(nrow(k$results) > 0 && all(k$results$to %in% c("C", "D", "E")))
  s <- spicy_glm(effect = "ratio", cells, condition = "response", r = 50, sigma = 50, test = "moderated")
  expect_true(nrow(s$results) > 0 && all(is.finite(s$results$p_value)))
})
