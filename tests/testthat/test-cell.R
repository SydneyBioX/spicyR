## spicyR Cell (method = "cell"): the C++ core against the plain-R reference written from the Supplementary's
## definitions (helper-reference.R), and the package's interface.

cells <- sim_cells()
types <- sort(unique(cells$cellType))
grid <- expand.grid(from = types, to = types, stringsAsFactors = FALSE)

test_that("the cell method matches the plain-R reference for every pair", {
  for (lc in c(FALSE, TRUE)) for (vr in c("cr2", "hartung_knapp")) {
    res <- suppressMessages(spicy(cells, "condition", subject = "patient", r = 30, availability = FALSE,
                                  labelClustering = lc, variance = vr))$cellResults
    # the reference is in (counted, centre) order: the user's from -> to is ref (to, from)
    ref <- ref_spicy_cell(transform(cells, subject = patient), data.frame(from = grid$to, to = grid$from), 30, lc, vr)
    ref <- data.frame(from = ref$to, to = ref$from, p_ref = ref$p_value, d_ref = ref$excess_difference)
    m <- merge(res, ref, by = c("from", "to"))
    expect_equal(nrow(m), nrow(grid))
    expect_lt(max(abs(m$p_value - m$p_ref)), 1e-8)
    expect_lt(max(abs(m$excess_difference - m$d_ref)), 1e-8)
  }
})

test_that("from is the centre and to the counted type", {
  res <- suppressMessages(spicy(cells, "condition", subject = "patient", r = 15, from = "tumour", to = "T"))$cellResults
  # T cells were moved next to tumour cells in group B: more T around each tumour cell
  expect_gt(res["tumour__T", "excess_difference"], 0)
  expect_lt(res["tumour__T", "p_value"], 0.01)
})

test_that("random labelling moments are exact (enumeration)", {
  set.seed(3)
  z <- data.frame(x = runif(14, 0, 60), y = runif(14, 0, 60), cellType = sample(c("A", "B", "C"), 14, TRUE))
  nb <- ref_neighbours(z$x, z$y, 25); ic <- ref_image_counts(z$cellType, nb, "A", "B")
  en <- ref_rl_enumerate(ic$ck, ic$n)
  expect_equal(unname(en[["mean"]]), ic$E, tolerance = 1e-12)
  expect_equal(unname(en[["var"]]), ic$v, tolerance = 1e-12)
})

ctx <- spicyR:::.cell_context(cells, "condition", "patient", "imageID", "cellType", c("x", "y"))
g <- spicyR:::.cell_graph(ctx, list(c("T", "tumour"), c("B", "macro")), r = 30)
rows <- spicyR:::.cell_rows(ctx, g, "T", "tumour")
m_units <- length(ctx$unit_labels)

test_that("the design path with the two group indicators equals the closed form", {
  a <- spicyR:::stats_excess_test(rows, rows$unit, rows$group, m_units, TRUE, "cr2")
  Z <- cbind(rows$group == 0, rows$group == 1)
  d <- spicyR:::stats_design_test(rows, rows$unit, m_units, Z, c(-1, 1), -1)
  expect_equal(d$tau2, a$tau2, tolerance = 1e-9)
  expect_equal(d$estimate, a$difference, tolerance = 1e-9)
  expect_equal(d$se, a$se, tolerance = 1e-8)
  expect_equal(d$df, a$df, tolerance = 1e-8)
})

test_that("covariates: the design test matches the dense definitions", {
  set.seed(4); zp <- rnorm(m_units)[rows$unit + 1L]; zi <- rnorm(nrow(rows))
  for (z in list(zp, zi)) {
    Z <- cbind(rows$group == 0, rows$group == 1, z - mean(z))
    a <- spicyR:::stats_design_test(rows, rows$unit, m_units, Z, c(-1, 1, 0), 0.01)
    b <- ref_design(rows, rows$unit, Z, c(-1, 1, 0), 0.01)
    expect_equal(c(a$estimate, a$se, a$df, a$p), unname(b), tolerance = 1e-8)
  }
})

test_that("Cox (Efron ties) matches survival::coxph", {
  set.seed(5); n <- 80; X <- cbind(rnorm(n), rbinom(n, 1, 0.4))
  tm <- round(stats::rexp(n, exp(0.4 * X[, 1])) * 5) / 5 + 0.2; ev <- rbinom(n, 1, 0.7)
  f <- survival::coxph(survival::Surv(tm, ev) ~ X); cc <- spicyR:::stats_cox_fit(tm, ev, X)
  expect_equal(cc$beta, unname(stats::coef(f)), tolerance = 1e-10)
  expect_equal(cc$se, unname(sqrt(diag(stats::vcov(f)))), tolerance = 1e-10)
  expect_equal(cc$martingale, unname(stats::residuals(f, "martingale")), tolerance = 1e-10)
  c0 <- spicyR:::stats_cox_fit(tm, ev, matrix(0, n, 0))
  expect_equal(c0$martingale, unname(stats::residuals(survival::coxph(survival::Surv(tm, ev) ~ 1), "martingale")),
               tolerance = 1e-10)
})

test_that("the survival score test is the design test on the martingale residuals", {
  set.seed(6); tm <- stats::rexp(m_units); ev <- rbinom(m_units, 1, 0.8)
  M <- spicyR:::stats_cox_fit(tm, ev, matrix(0, m_units, 0))$martingale
  s <- spicyR:::stats_survival_test(rows, rows$unit, m_units, M, tm, ev)
  d <- spicyR:::stats_design_test(rows, rows$unit, m_units, cbind(1, M[rows$unit + 1L]), c(0, 1), -1)
  expect_true(s$ok)
  expect_equal(s$score_p, d$p, tolerance = 1e-12)
  pheno <- unique(cells[, c("patient", "imageID")])
  cells2 <- cells; i <- match(cells2$patient, unique(cells$patient))
  cells2$os <- survival::Surv(tm[i], ev[i])
  res <- suppressMessages(spicy(cells2, "os", subject = "patient", r = 30, from = "tumour", to = c("T", "B")))
  expect_s4_class(res, "SpicyResults")
  expect_true(all(c("hazard_ratio_sd", "p_value") %in% names(res$cellResults)))
})

test_that("max-T: one radius gives the single-radius p, and the tail matches an exact integral", {
  expect_equal(spicyR:::stats_max_t(matrix(rnorm(10), 10), 2.1, 7)$p, 2 * stats::pt(-2.1, 7), tolerance = 1e-12)
  r <- 0.8; cv <- 3.5
  f <- function(x) stats::dnorm(x) * (stats::pnorm((-cv - r * x) / sqrt(1 - r^2)) + stats::pnorm((-cv + r * x) / sqrt(1 - r^2)))
  exact <- 2 * stats::pnorm(-cv) + stats::integrate(f, -cv, cv, rel.tol = 1e-12)$value
  expect_equal(spicyR:::stats_mvn_outside(c(-cv, -cv), c(cv, cv), matrix(c(1, r, r, 1), 2)), exact, tolerance = 1e-8)
  res <- suppressMessages(spicy(cells, "condition", subject = "patient", r = c(10, 20, 40), from = "tumour", to = "T"))
  expect_true(is.finite(res$cellResults$p_value))
  expect_equal(nrow(res$radiusResults), 3L)
})

test_that("results work with topPairs, bind, spicyBoxPlot and signifPlot", {
  res <- suppressMessages(spicy(cells, "condition", subject = "patient", r = 30))
  tp <- topPairs(res, n = 5)
  expect_equal(nrow(tp), 5L)
  expect_true(all(c("from", "to", "coefficient", "p.value") %in% names(tp)))
  expect_equal(ncol(bind(res)), 3L + nrow(res$cellResults))
  expect_s3_class(spicyBoxPlot(res, from = "tumour", to = "T"), "ggplot")
  expect_no_error(signifPlot(res))
})

test_that("more than two conditions give one contrast per level", {
  c3 <- cells; c3$condition <- as.character(c3$condition)
  c3$condition[c3$patient %in% c("A01", "A02", "A03")] <- "C"
  res <- suppressMessages(spicy(c3, "condition", subject = "patient", r = 30, from = "tumour", to = "T"))
  expect_equal(sort(unique(res$cellResults$level)), c("B", "C"))
  expect_equal(colnames(res$p.value), c("(Intercept)", "conditionB", "conditionC"))
})
