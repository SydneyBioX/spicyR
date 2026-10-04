## spicyR Cell (method = "cell"): the C++ core against the plain-R reference written from the Supplementary's
## definitions (helper-reference.R), and the package's interface.

cells <- sim_cells()
types <- sort(unique(cells$cellType))
grid <- expand.grid(from = types, to = types, stringsAsFactors = FALSE)

test_that("the cell method matches the plain-R reference for every pair", {
  for (ef in c("allocation", "count")) for (lc in c(FALSE, TRUE)) for (vr in c("cr2", "hartung_knapp")) {
    res <- suppressMessages(spicy(cells, "condition", subject = "patient", r = 30, effect = ef,
                                  labelClustering = lc, variance = vr))$cellResults
    # the reference is in (counted, centre) order, which is the user's from -> to
    ref <- ref_spicy_cell(transform(cells, subject = patient), grid, 30, lc, vr, effect = ef)
    ref <- data.frame(from = ref$from, to = ref$to, p_ref = ref$p_value, d_ref = ref$excess_difference)
    m <- merge(res, ref, by = c("from", "to"))
    expect_equal(nrow(m), nrow(grid))
    expect_lt(max(abs(m$p_value - m$p_ref)), 1e-8)
    expect_lt(max(abs(m$excess_difference - m$d_ref)), 1e-8)
  }
})

test_that("to is the centre (relabelled) type and from the counted type", {
  for (ef in c("allocation", "count")) {
    res <- suppressMessages(spicy(cells, "condition", subject = "patient", r = 15, from = "tumour", to = "T",
                                  effect = ef))$cellResults
    # T cells were moved next to tumour cells in group B: more T cells have a tumour cell nearby, and there are more
    # tumour cells around each T cell
    expect_gt(res["tumour__T", "excess_difference"], 0)
    expect_lt(res["tumour__T", "p_value"], 0.01)
  }
})

test_that("allocation: the effect is the fraction of T cells moved next to tumour", {
  # 40% of group B's T cells were moved within 6 units of a tumour cell; at r = 15 nearly all of them have a tumour cell
  # within r, so the extra fraction is a little under 0.4 (the moved cells that already had a tumour cell nearby are
  # not extra)
  res <- suppressMessages(spicy(cells, "condition", subject = "patient", r = 15, from = "tumour", to = "T"))
  expect_identical(res$effect, "allocation")
  x <- res$cellResults["tumour__T", ]
  expect_lt(abs(x$excess_ref), 0.05)
  expect_gt(x$excess_comp, 0.25); expect_lt(x$excess_comp, 0.45)
})

test_that("allocation moments are exact: cross-pairs by enumeration, self-pairs by enumeration of every labelling", {
  set.seed(3)
  z <- data.frame(x = runif(14, 0, 60), y = runif(14, 0, 60), cellType = sample(c("A", "B", "C"), 14, TRUE))
  nb <- ref_neighbours(z$x, z$y, 25)
  ic <- ref_image_counts_alloc(z$cellType, nb, "A", "B")
  en <- ref_rl_enumerate(ic$ck, ic$n_to)
  expect_equal(unname(en[["mean"]]), ic$E, tolerance = 1e-12)
  expect_equal(unname(en[["var"]]), ic$v, tolerance = 1e-12)
  # self-pair: every placement of the A labels among the 14 cells
  ia <- ref_image_counts_alloc(z$cellType, nb, "A", "A"); nA <- sum(z$cellType == "A")
  S <- utils::combn(14, nA)
  O <- apply(S, 2, function(s) sum(rowSums(nb$adj[s, s, drop = FALSE]) > 0))
  expect_equal(mean(O), ia$E, tolerance = 1e-12)
  # and the core's expectation
  ctx1 <- spicyR:::.cell_context(transform(z, imageID = "i", condition = "a"), NULL, NULL, "imageID", "cellType",
                                 c("x", "y"))
  spicyR:::dataset_build_radius_index(ctx1$data, 25)
  se <- spicyR:::dataset_self_any_expected(ctx1$data, FALSE, length(ctx1$type_labels))
  expect_equal(se[match("A", ctx1$type_labels), 1], ia$E, tolerance = 1e-10)
})

test_that("k nearest neighbours: both effects run, and images with at most k cells are skipped", {
  ck <- rbind(cells, data.frame(x = c(1, 2, 3), y = c(1, 2, 3), cellType = c("tumour", "T", "B"),
                                imageID = "tiny", patient = "A01", condition = factor("A", levels = c("A", "B"))))
  for (ef in c("allocation", "count")) {
    res <- suppressMessages(spicy(ck, "condition", subject = "patient", k = 5, from = "tumour", to = "T",
                                  effect = ef))$cellResults
    expect_gt(res["tumour__T", "excess_difference"], 0)
    expect_lt(res["tumour__T", "p_value"], 0.05)
  }
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
  s <- spicyR:::stats_survival_test(rows, rows$unit, m_units, M, tm, ev, numeric(0))
  d <- spicyR:::stats_design_test(rows, rows$unit, m_units, cbind(1, M[rows$unit + 1L]), c(0, 1), -1)
  expect_true(s$ok)
  expect_equal(s$score_p, d$p, tolerance = 1e-12)
  # adjusted for abundance: the log share added to the design, tau2 held at the unadjusted value
  x <- spicyR:::.cell_share(ctx, rows, "T")
  sa <- spicyR:::stats_survival_test(rows, rows$unit, m_units, M, tm, ev, x)
  da <- spicyR:::stats_design_test(rows, rows$unit, m_units, cbind(1, M[rows$unit + 1L], x - mean(x)), c(0, 1, 0), d$tau2)
  expect_equal(sa$score_p, da$p, tolerance = 1e-12)
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

test_that("adjustAbundance = TRUE adjusts for abundance (the dense definition), with the unadjusted test alongside", {
  # rows are the core's (counted T, centre tumour), i.e. the pair T -> tumour; the covariate is the share of T
  share <- spicyR:::.cell_share(ctx, rows, "T")
  Z <- cbind(rows$group == 0, rows$group == 1, share - mean(share))
  a <- spicyR:::stats_excess_test(rows, rows$unit, rows$group, m_units, TRUE, "cr2")
  b <- ref_design(rows, rows$unit, Z, c(-1, 1, 0), a$tau2)
  res <- suppressMessages(spicy(cells, "condition", subject = "patient", r = 30, from = "T", to = "tumour",
                                adjustAbundance = TRUE))$cellResults
  expect_equal(c(res$excess_difference, res$se, res$df, res$p_value), unname(b), tolerance = 1e-8)
  expect_equal(res$adjusted_for, "abundance")
  expect_equal(res$unadjusted_p_value, a$p, tolerance = 1e-10)
  # the default is unadjusted
  res0 <- suppressMessages(spicy(cells, "condition", subject = "patient", r = 30, from = "T", to = "tumour"))$cellResults
  expect_equal(res0$p_value, res$unadjusted_p_value, tolerance = 1e-12)
  expect_null(res0$unadjusted_p_value)
})

test_that("several contrasts of one design: influences give the variance, and each covariate has its effect", {
  set.seed(5); zp <- rnorm(m_units)[rows$unit + 1L]; zi <- rnorm(nrow(rows))
  Z <- cbind(rows$group == 0, rows$group == 1, zp - mean(zp), zi - mean(zi))
  C <- rbind(c(-1, 1, 0, 0), c(0, 0, 1, 0), c(0, 0, 0, 1))
  d <- spicyR:::stats_design_tests(rows, rows$unit, m_units, Z, C, 0.01, FALSE)
  for (j in 1:3) {
    one <- spicyR:::stats_design_test(rows, rows$unit, m_units, Z, C[j, ], 0.01)
    expect_equal(c(d[[j]]$estimate, d[[j]]$se, d[[j]]$df), c(one$estimate, one$se, one$df), tolerance = 1e-10)
    expect_equal(sum(d[[j]]$influence^2), d[[j]]$se^2, tolerance = 1e-10)
  }
  cz <- cells; cz$age <- zp[match(cz$patient, ctx$unit_labels[rows$unit + 1L])]
  cz$batch <- c("a", "b", "c")[(match(cz$imageID, unique(cz$imageID)) %% 3) + 1]
  res <- suppressMessages(spicy(cz, "condition", subject = "patient", r = 30, from = "tumour", to = "T",
                                covariates = c("age", "batch")))$cellResults
  expect_true(all(c("age_effect", "age_p_value", "batchb_effect", "batchc_p_value") %in% names(res)))
  expect_equal(res$adjusted_for, "covariates")
  res <- suppressMessages(spicy(cz, "condition", subject = "patient", r = 30, from = "tumour", to = "T",
                                covariates = c("age", "batch"), adjustAbundance = TRUE))$cellResults
  expect_true("abundance_effect" %in% names(res))
  expect_equal(res$adjusted_for, "abundance+covariates")
})

test_that("image weights sum to one within each condition", {
  res <- suppressMessages(spicy(cells, "condition", subject = "patient", r = 30, from = "tumour", to = "T"))
  w <- res$imageWeights[["tumour__T"]]
  expect_equal(as.vector(tapply(w, res$condition, sum, na.rm = TRUE)), c(1, 1), tolerance = 1e-10)
})

test_that("a pair's result does not depend on the other pairs requested (psi over every counted type)", {
  for (ef in c("allocation", "count")) {
    full <- suppressMessages(spicy(cells, "condition", subject = "patient", r = 30, effect = ef,
                                   adjustAbundance = TRUE))$cellResults
    one <- suppressMessages(spicy(cells, "condition", subject = "patient", r = 30, from = "tumour", to = "T", effect = ef,
                                  adjustAbundance = TRUE))$cellResults
    for (col in c("p_value", "excess_difference", "unadjusted_p_value"))
      expect_equal(one["tumour__T", col], full["tumour__T", col], tolerance = 1e-12)
  }
})

test_that("allocation, avoidance: the effect is minus the fraction moved away, however many from cells there are", {
  # Group A has few `from` cells, group B many. In both, 30% of the `to` cells that random labelling would put next to
  # a `from` cell are placed away from them instead, so the effect is -0.3 in both groups and the difference is ~0.
  set.seed(11)
  one_image <- function(id, patient, cond, n_from, f = 0.3, N = 900, n_to = 150, r = 4) {
    xy <- cbind(runif(N, 0, 100), runif(N, 0, 100))
    type <- rep("other", N); type[sample(N, n_from)] <- "A"
    cand <- which(type != "A")
    D <- as.matrix(stats::dist(xy))
    near <- cand[vapply(cand, function(i) any(D[i, type == "A"] <= r), TRUE)]
    far <- setdiff(cand, near)
    q <- length(near) / length(cand)
    k_near <- round(n_to * q * (1 - f))
    type[c(sample(near, k_near), sample(far, n_to - k_near))] <- "T"
    data.frame(x = xy[, 1], y = xy[, 2], cellType = type, imageID = id, patient = patient,
               condition = factor(cond, levels = c("A", "B")))
  }
  ck <- do.call(rbind, c(lapply(1:4, function(i) one_image(paste0("a", i), paste0("pa", i), "A", 25)),
                         lapply(1:4, function(i) one_image(paste0("b", i), paste0("pb", i), "B", 120))))
  x <- suppressMessages(spicy(ck, "condition", subject = "patient", r = 4, from = "A", to = "T",
                              labelClustering = FALSE))$cellResults["A__T", ]
  expect_identical(x$side, "avoid")
  expect_lt(abs(x$excess_ref + 0.3), 0.03)
  expect_lt(abs(x$excess_comp + 0.3), 0.03)
  expect_gt(x$p_value, 0.05)
})

test_that("variance = 'auto' uses Hartung-Knapp when a condition has at most 5 patients, CR2 otherwise", {
  pts <- unique(cells[, c("patient", "condition")])
  m_min <- min(table(pts$condition))
  res <- suppressMessages(spicy(cells, "condition", subject = "patient", r = 15, from = "tumour", to = "T"))
  expect_identical(res$variance, if (m_min <= 5) "hartung_knapp" else "cr2")
  small <- cells[cells$patient %in% unlist(lapply(split(pts$patient, pts$condition), head, 4)), ]
  rs <- suppressMessages(spicy(small, "condition", subject = "patient", r = 15, from = "tumour", to = "T"))
  expect_identical(rs$variance, "hartung_knapp")
  hk <- suppressMessages(spicy(small, "condition", subject = "patient", r = 15, from = "tumour", to = "T",
                               variance = "hartung_knapp"))
  expect_equal(rs$cellResults$p_value, hk$cellResults$p_value)
})
