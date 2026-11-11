## Survival outcomes with the cell method (DECISIONS 5 Nov 2026 and 7 Nov 2026): patients with no outcome left out
## first, covariates on both sides, one radius, the "enough information" rules, and the df guard of every test.

sv_cells <- function() {
  d <- sim_cells(seed = 11, n_pat = 7)   # 14 patients, two images each
  pats <- sort(unique(d$patient))
  set.seed(12)
  tm <- stats::setNames(round(stats::rexp(length(pats), 0.1), 1) + 0.1, pats)
  ev <- stats::setNames(c(1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1), pats)
  d$os <- survival::Surv(tm[d$patient], ev[d$patient])
  d$age <- stats::setNames(round(stats::rnorm(length(pats), 60, 10)), pats)[d$patient]           # patient-level
  d$depth <- stats::setNames(stats::rnorm(length(unique(d$imageID))), unique(d$imageID))[d$imageID]  # image-level
  d
}
scells <- sv_cells()
spats <- sort(unique(scells$patient))

## One pair's rows through the package's own steps (no condition: every image), and the unit codes' outcome.
sv_rows <- function(d, f = "tumour", t = "T", r = 30) {
  ctx <- spicyR:::.cell_context(d, NULL, "patient", "imageID", "cellType", c("x", "y"), survival = TRUE)
  g <- spicyR:::.cell_graph(ctx, list(c(f, t)), r = r, label_clustering = FALSE)
  g$image_heterogeneity <- TRUE
  first <- match(ctx$unit_labels, d$patient)
  list(ctx = ctx, rows = spicyR:::.cell_rows(ctx, g, f, t), time = d$os[first, 1], event = as.integer(d$os[first, 2]),
       age = d$age[first], m = length(ctx$unit_labels))
}

test_that("every pair has a row; hr_p_value is gone; n_deaths and not_tested are reported", {
  res <- spicy(scells, condition = "os", subject = "patient", r = 30, from = "tumour")$cellResults
  expect_equal(nrow(res), length(unique(scells$cellType)))
  expect_false("hr_p_value" %in% names(res))
  expect_true(all(c("n_deaths", "not_tested", "log_hr_se", "hazard_ratio_sd") %in% names(res)))
  tested <- is.na(res$not_tested)
  expect_true(any(tested))
  expect_true(all(is.finite(res$p_value[tested])))
  expect_true(all(is.na(res$p_value[!tested])))
  expect_true(all(res$n_deaths <= sum(c(1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1))))
})

test_that("a survival outcome takes one radius: several are an error", {
  expect_error(spicy(scells, condition = "os", subject = "patient", r = c(20, 30), from = "tumour", to = "T"),
               "one radius")
  expect_error(spicy(scells, condition = "os", subject = "patient", r = c(10, 20, 30, 40, 50, 60, 70), from = "tumour",
                     to = "T"), "one radius")
})

test_that("patients with no outcome are left out first: the result is that of the data without them", {
  d <- scells; gone <- spats[c(2, 5, 9)]
  tm <- d$os[, 1]; tm[d$patient %in% gone] <- NA
  d$os <- survival::Surv(tm, d$os[, 2])
  without <- scells[!scells$patient %in% gone, ]
  for (ef in c("allocation", "count")) for (adj in c(FALSE, TRUE)) for (cv in list(NULL, c("age", "depth"))) {
    a <- spicy(d, condition = "os", subject = "patient", r = 30, from = c("tumour", "T"), effect = ef,
               adjustAbundance = adj, covariates = cv)$cellResults
    b <- spicy(without, condition = "os", subject = "patient", r = 30, from = c("tumour", "T"), effect = ef,
               adjustAbundance = adj, covariates = cv)$cellResults
    expect_identical(a$not_tested, b$not_tested)
    num <- names(a)[vapply(a, is.numeric, TRUE)]
    expect_equal(a[num], b[num], tolerance = 1e-8)
    if (ef == "allocation") expect_identical(a$side, b$side)
  }
})

test_that("covariates are on both sides: in the null Cox model and in the fraction model", {
  s <- sv_rows(scells, "T", "tumour"); rows <- s$rows; m <- s$m   # a pair with tau2 > 0, so it has a hazard ratio
  W <- cbind(age = s$age)
  M <- spicyR:::stats_cox_fit(s$time, s$event, W)$martingale
  expect_equal(M, unname(stats::residuals(survival::coxph(survival::Surv(s$time, s$event) ~ W), "martingale")),
               tolerance = 1e-8)
  age_c <- s$age[rows$unit + 1L] - mean(s$age[rows$unit + 1L])
  z <- spicyR:::stats_survival_test(rows, rows$unit, m, M, s$time, s$event, cbind(age = s$age[rows$unit + 1L]), W)
  expect_true(z$ok)
  # the score test: the M coefficient of (1, M, age), tau2 by Paule-Mandel under that (patient-level) design
  d <- spicyR:::stats_design_test(rows, rows$unit, m, cbind(1, M[rows$unit + 1L], age_c), c(0, 1, 0), -1)
  expect_equal(z$score_p, d$p, tolerance = 1e-12)
  expect_equal(z$score_df, d$df, tolerance = 1e-12)
  # the hazard ratio: the shrunken excess of (1, age) (without M), at the average age, in a Cox model with age
  h <- spicyR:::stats_design_test(rows, rows$unit, m, cbind(1, age_c), c(1, 0), -1)
  y <- rows$O - rows$E - rows$n * h$theta[2] * age_c
  sk <- tapply(rows$n * y / rows$v, rows$unit, sum); J <- tapply(rows$n^2 / rows$v, rows$unit, sum)
  blup <- h$theta[1] + h$tau2 / (h$tau2 + 1 / J) * (sk / J - h$theta[1])
  u <- as.integer(names(J)) + 1L
  hs <- (blup - mean(blup)) / stats::sd(blup)
  cx <- survival::coxph(survival::Surv(s$time[u], s$event[u]) ~ hs + s$age[u])
  expect_equal(z$log_hr_sd, unname(stats::coef(cx)[1]), tolerance = 1e-6)
  expect_equal(z$hr_se, unname(sqrt(stats::vcov(cx)[1, 1])), tolerance = 1e-6)
  expect_equal(z$tau2, h$tau2, tolerance = 1e-12)
  expect_gt(z$tau2, 0)
  # an image-level covariate: tau2 held at that of the patient-level columns (1, M); TODO(ext-AB): option (iv)
  dep <- scells$depth[match(s$ctx$image_labels[rows$img + 1L], scells$imageID)]
  z2 <- spicyR:::stats_survival_test(rows, rows$unit, m, M, s$time, s$event, cbind(dep), matrix(0, m, 0))
  d0 <- spicyR:::stats_design_test(rows, rows$unit, m, cbind(1, M[rows$unit + 1L]), c(0, 1), -1)
  d2 <- spicyR:::stats_design_test(rows, rows$unit, m, cbind(1, M[rows$unit + 1L], dep - mean(dep)), c(0, 1, 0), d0$tau2)
  expect_equal(z2$score_p, d2$p, tolerance = 1e-12)
  expect_identical(z2$q, 4L)   # 1, M, depth, and the patient mean depth will have
  # through spicy(): the same test
  res <- spicy(scells, condition = "os", subject = "patient", r = 30, from = "T", to = "tumour",
               covariates = "age")$cellResults
  expect_equal(res$p_value, z$score_p, tolerance = 1e-10)
  expect_equal(res$hazard_ratio_sd, exp(z$log_hr_sd), tolerance = 1e-10)
  expect_identical(res$adjusted_for, "covariates")
})

test_that("enough information: a death and a survivor, and at least q + 2 informative patients", {
  s <- sv_rows(scells); rows <- s$rows; m <- s$m
  none <- function(n) matrix(0, n, 0)
  test <- function(keep, event = s$event, X = none(nrow(rows)), W = none(m)) {
    M <- rep(NA_real_, m)
    M[keep] <- spicyR:::stats_cox_fit(s$time[keep], event[keep], W[keep, , drop = FALSE])$martingale
    spicyR:::stats_survival_test(rows, rows$unit, m, M, s$time, event, X, W)
  }
  all_units <- seq_len(m)
  z <- test(all_units, event = rep(0L, m))
  expect_false(z$ok); expect_identical(z$reason, "no_deaths"); expect_identical(z$n_deaths, 0L)
  z <- test(all_units, event = rep(1L, m))
  expect_false(z$ok); expect_identical(z$reason, "no_survivors")
  # without covariates q = 2: 4 informative patients are tested, 3 are not
  four <- c(which(s$event == 1)[1:2], which(s$event == 0)[1:2])
  z <- test(four[1:3])
  expect_identical(z$reason, "too_few_units"); expect_identical(z$n_units, 3L)
  z <- test(four)
  expect_identical(z$n_units, 4L); expect_identical(z$q, 2L); expect_identical(z$n_deaths, 2L)
  expect_false(identical(z$reason, "too_few_units"))
  # one patient-level covariate (q = 3): 5 are needed
  W <- cbind(age = s$age); X <- cbind(age = s$age[rows$unit + 1L])
  expect_identical(test(four, X = X, W = W)$reason, "too_few_units")
  five <- c(four, which(s$event == 1)[3])
  expect_false(identical(test(five, X = X, W = W)$reason, "too_few_units"))
  # through spicy(): a pair with no deaths among its patients
  d <- scells; d$os <- survival::Surv(d$os[, 1], rep(0, nrow(d)))
  res <- spicy(d, condition = "os", subject = "patient", r = 30, from = "tumour", to = c("T", "B"))$cellResults
  expect_identical(res$not_tested, c("no_deaths", "no_deaths"))
  expect_true(all(is.na(res$p_value)))
})

test_that("the df guard: a test on fewer than 1 Satterthwaite df, or infinite df, is not done", {
  rules <- spicyR:::stats_test_rules(TRUE)
  expect_identical(rules$min_df, 1)
  expect_true(spicyR:::.cell_too_few_df(0.7, rules))
  expect_true(spicyR:::.cell_too_few_df(Inf, rules))
  expect_true(spicyR:::.cell_too_few_df(NaN, rules))
  expect_false(spicyR:::.cell_too_few_df(1, rules))
  # four patients, one with nearly all the weight and the one death: the score test's df are 0.03
  rows <- data.frame(img = 0:3, O = c(-0.986191, -0.394855, 1.30852, 0.95476), E = 0, n = 1,
                     v = c(0.000196998, 71.9433, 0.00535398, 8.92664))
  M <- c(4.89172, -0.127562, -0.0564688, -0.051067)
  d <- spicyR:::stats_design_test(rows, 0:3, 4L, cbind(1, M), c(0, 1), -1)
  expect_true(d$ok); expect_lt(d$df, 1)
  z <- spicyR:::stats_survival_test(rows, 0:3, 4L, M, 1:4, c(1L, 0L, 0L, 0L), matrix(0, 4, 0), matrix(0, 4, 0))
  expect_false(z$ok); expect_identical(z$reason, "too_few_df")
})
