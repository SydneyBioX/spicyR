## spicyR Cell (method = "cell"): the C++ core against the plain-R reference written from the Supplementary's
## definitions (helper-reference.R), and the package's interface.

cells <- sim_cells()
types <- sort(unique(cells$cellType))
grid <- expand.grid(from = types, to = types, stringsAsFactors = FALSE)

## The cell method through its internal steps, with the label-clustering inflation available (spicy() has it off)
cell_internal <- function(cells, r, effect, label_clustering, variance, image_heterogeneity = TRUE) {
  ctx <- spicyR:::.cell_context(cells, "condition", "patient", "imageID", "cellType", c("x", "y"), dataset = TRUE)
  pairs <- lapply(seq_len(nrow(grid)), function(i) c(grid$from[i], grid$to[i]))
  g <- spicyR:::.cell_graph(ctx, pairs, r = r, label_clustering = label_clustering, effect = effect)
  g$image_heterogeneity <- image_heterogeneity
  fits <- lapply(pairs, function(p) spicyR:::.cell_pair_test(ctx, g, p[1], p[2], TRUE, variance, FALSE))
  spicyR:::.cell_table(fits, ctx, FALSE)
}

test_that("spicy() is the internal cell method with the frailty and image random effect, without inflation", {
  for (ef in c("allocation", "count")) {
    res <- suppressMessages(spicy(cells, "condition", subject = "patient", r = 30, effect = ef))$cellResults
    int <- cell_internal(cells, 30, ef, FALSE, "cr2")
    m <- merge(res, int, by = c("from", "to"))
    expect_equal(nrow(m), nrow(grid))
    expect_equal(m$p_value.x, m$p_value.y, tolerance = 1e-12)
  }
  expect_error(spicy(cells, "condition", subject = "patient", labelClustering = TRUE), "not an argument")
})

test_that("the cell method matches the plain-R reference for every pair", {
  for (ef in c("allocation", "count")) for (lc in c(FALSE, TRUE)) for (vr in c("cr2", "hartung_knapp")) {
    res <- cell_internal(cells, 30, ef, lc, vr)
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
                                 c("x", "y"), dataset = TRUE)
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
g <- spicyR:::.cell_graph(ctx, list(c("T", "tumour"), c("B", "macro")), r = 30, label_clustering = FALSE)
rows_raw <- spicyR:::.cell_rows(ctx, g, "T", "tumour")
g$image_heterogeneity <- TRUE   # as spicy() by default
rows <- spicyR:::.cell_rows(ctx, g, "T", "tumour")

test_that("images within a patient: the core's rows match the reference", {
  ref <- ref_image_heterogeneity(rows_raw, rows_raw$unit)
  expect_equal(attr(rows, "tau2_image"), attr(ref, "tau2_image"), tolerance = 1e-12)
  for (k in c("O", "E", "n", "v")) expect_equal(rows[[k]], ref[[k]], tolerance = 1e-12)
  # the patient's images agree here (no planted spread), and two images per patient carry a positive tau2 only
  # when their effects differ beyond the sampling variance: check one where they do
  r2 <- rows_raw; r2$O <- r2$O + ifelse(r2$img %% 2 == 0, 1, -1) * 3 * sqrt(r2$v)
  expect_gt(attr(spicyR:::stats_image_heterogeneity_rows(r2, r2$unit), "tau2_image"), 0)
  expect_equal(attr(spicyR:::stats_image_heterogeneity_rows(r2, r2$unit), "tau2_image"),
               attr(ref_image_heterogeneity(r2, r2$unit), "tau2_image"), tolerance = 1e-12)
})

test_that("with one image per patient, the image random effect changes nothing", {
  one <- cells[grepl("_1$", cells$imageID), ]
  a <- cell_internal(one, 30, "allocation", FALSE, "cr2", TRUE)
  b <- cell_internal(one, 30, "allocation", FALSE, "cr2", FALSE)
  expect_equal(a$p_value, b$p_value, tolerance = 1e-10)
  expect_equal(a$excess_difference, b$excess_difference, tolerance = 1e-10)
  expect_true(all(a$tau2_image == 0))
})
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

test_that("max-T: the Genz-Bretz quasi-Monte Carlo matches mvtnorm::pmvnorm and repeats exactly", {
  skip_if_not_installed("mvtnorm")
  ar1 <- function(K, rho) rho^abs(outer(seq_len(K), seq_len(K), "-"))
  mats <- list(ar1(3, 0.6), ar1(4, 0.9), ar1(6, 0.5),
               stats::cov2cor(crossprod(matrix(c(3, 1, 0, 2, 1, 2, 2, 0, 1, 1, 3, 1, 0, 2, 1, 2, 1, 1, 0, 3), 5))))
  cvs <- c(2.2, 2.8, 3.1, 2.5)
  for (i in seq_along(mats)) {
    R <- mats[[i]]; K <- nrow(R); cv <- rep(cvs[i], K)
    p <- spicyR:::stats_mvn_outside(-cv, cv, R)
    set.seed(1)
    ref <- 1 - mvtnorm::pmvnorm(-cv, cv, corr = R, algorithm = mvtnorm::GenzBretz(maxpts = 2e6, abseps = 1e-6, releps = 0))[1]
    expect_lt(abs(p - ref), 2e-4)
    expect_identical(spicyR:::stats_mvn_outside(-cv, cv, R), p)
    q <- spicyR:::stats_mvn_outside_qmc(-cv, cv, R)
    expect_identical(q$p, p)
    expect_lt(q$error, 1e-4)
  }
})

test_that("Davies: one radius gives the single-radius p, matches a plain-R reference and bounds the max-T", {
  expect_equal(spicyR:::stats_davies(matrix(rnorm(10), 10), 2.1, 7)$p, 2 * stats::pt(-2.1, 7), tolerance = 1e-12)
  set.seed(3); m <- 40; K <- 5
  # nested radii: each radius adds neighbours to the previous one, so adjacent radii are highly correlated
  base <- rnorm(m); U <- t(apply(matrix(rnorm(m * K), m), 1, cumsum)) + 2 * base
  tt <- c(1.2, 2.4, 2.9, 2.6, 1.8); df <- rep(30, K)
  R <- stats::cov2cor(crossprod(U)); z <- stats::qnorm(stats::pt(-abs(tt), df), lower.tail = FALSE); cc <- max(z)
  V <- sum(acos(R[cbind(1:(K - 1), 2:K)]))
  ref <- min(1, 2 * stats::pnorm(-cc) + V * exp(-cc^2 / 2) / pi)
  d <- spicyR:::stats_davies(U, tt, df)
  expect_equal(d$p, ref, tolerance = 1e-10)
  expect_equal(d$best, 3L)
  mt <- spicyR:::stats_max_t(U, tt, df)
  expect_gte(d$p, mt$p * (1 - 1e-6))
  expect_lt(d$p / mt$p, 1.25)
  res <- suppressMessages(spicy(cells, "condition", subject = "patient", r = seq(5, 40, by = 5), combine = "davies",
                                from = "tumour", to = "T"))
  expect_true(is.finite(res$cellResults$p_value))
  expect_equal(nrow(res$radiusResults), 8L)
  expect_error(spicy(cells, "condition", subject = "patient", r = seq(5, 40, by = 5), from = "tumour", to = "T"),
               "up to 6 radii")
})

test_that("several radii from one neighbour search equal the single-radius graphs", {
  ctx <- spicyR:::.cell_context(cells, "condition", "patient", "imageID", "cellType", c("x", "y"))
  radii <- c(5, 12.5, 20, 40)
  for (eff in c("allocation", "count")) {
    G <- spicyR:::.cell_graphs_radii(ctx, radii, eff)
    for (i in seq_along(radii)) {
      g1 <- spicyR:::.cell_graph(ctx, list(c("tumour", "T")), r = radii[i], label_clustering = FALSE, effect = eff)
      for (nm in if (eff == "allocation") c("any", "self_expected") else c("totals", "sq")) expect_identical(G[[i]][[nm]], g1[[nm]])
    }
  }
  res <- suppressMessages(spicy(cells, "condition", subject = "patient", r = radii, combine = "davies", from = "tumour", to = c("T", "B")))
  one <- lapply(radii, function(r) suppressMessages(spicy(cells, "condition", subject = "patient", r = r, from = "tumour", to = c("T", "B")))$cellResults)
  for (i in seq_along(radii)) {
    rr <- res$radiusResults[res$radiusResults$r == radii[i], ]
    expect_equal(rr$p_value[match(one[[i]]$to, rr$to)], one[[i]]$p_value, tolerance = 1e-12)
  }
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
  x <- suppressMessages(spicy(ck, "condition", subject = "patient", r = 4, from = "A", to = "T"))$cellResults["A__T", ]
  expect_identical(x$side, "avoid")
  expect_lt(abs(x$excess_ref + 0.3), 0.03)
  expect_lt(abs(x$excess_comp + 0.3), 0.03)
  expect_gt(x$p_value, 0.05)
})

test_that("variance is CR2 by default, whatever the number of patients; Hartung-Knapp on request", {
  pts <- unique(cells[, c("patient", "condition")])
  small <- cells[cells$patient %in% unlist(lapply(split(pts$patient, pts$condition), head, 4)), ]
  def <- suppressMessages(spicy(small, "condition", subject = "patient", r = 15, from = "tumour", to = "T"))
  expect_identical(def$variance, "cr2")
  hk <- suppressMessages(spicy(small, "condition", subject = "patient", r = 15, from = "tumour", to = "T",
                               variance = "hartung_knapp"))
  expect_identical(hk$variance, "hartung_knapp")
  expect_equal(hk$cellResults$df, 6)
  expect_error(spicy(small, "condition", subject = "patient", variance = "auto"))
})

test_that("images with a missing condition are left out of the covariate design, as of the unadjusted test", {
    set.seed(2)
    n <- 20
    im <- data.frame(imageID = sprintf("i%02d", 1:n), cond = factor(rep(c("a", "b"), n / 2)), cv = rnorm(n))
    im$cond[1:2] <- NA
    cells <- do.call(rbind, lapply(1:n, function(i) data.frame(imageID = im$imageID[i],
        cellType = sample(c("A", "B"), 500, TRUE), x = runif(500, 0, 500), y = runif(500, 0, 500))))
    for (effect in c("count", "allocation")) {
        r <- spicy(cells, condition = "cond", covariates = "cv", imageInfo = im, cellType = "cellType", r = 50,
                   effect = effect)
        expect_true(all(r$cellResults$adjusted_for == "covariates"))
    }
    keep <- !is.na(im$cond)
    full <- spicy(cells, condition = "cond", covariates = "cv", imageInfo = im, cellType = "cellType", r = 50,
                  effect = "count")
    sub <- spicy(cells[cells$imageID %in% im$imageID[keep], ], condition = "cond", covariates = "cv",
                 imageInfo = im[keep, ], cellType = "cellType", r = 50, effect = "count")
    a <- full$cellResults[order(full$cellResults$from, full$cellResults$to), ]
    b <- sub$cellResults[order(sub$cellResults$from, sub$cellResults$to), ]
    expect_equal(a$p_value, b$p_value)
    expect_equal(a$excess_difference, b$excess_difference)
})

## tumour -> T with no variation at two radii: no T cell has a tumour cell within 2 (B cells do), and every T cell has
## one within 10 (the far cells do not); at r = 4 the fraction differs by condition.
nv_cells <- function(n_img = 16, seed = 3) {
  set.seed(seed)
  do.call(rbind, lapply(seq_len(n_img), function(i) {
    g <- expand.grid(x = seq(10, 90, 20), y = seq(10, 90, 20)); nt <- sample(12:25, 1)
    g <- g[sample(nrow(g), nt), ]; tx <- g$x; ty <- g$y
    near <- function(n, lo, hi) { k <- sample(nt, n, TRUE); a <- runif(n, 0, 2 * pi); d <- runif(n, lo, hi)
      cbind(tx[k] + d * cos(a), ty[k] + d * sin(a)) }
    grp <- if (i <= n_img / 2) "A" else "B"
    nT <- sample(30:50, 1); nc <- round(nT * if (grp == "A") runif(1, 0.2, 0.5) else runif(1, 0.5, 0.8))
    nB <- sample(20:40, 1); nF <- sample(30:50, 1)
    xy <- rbind(cbind(tx, ty), near(nc, 2.5, 3.5), near(nT - nc, 4.5, 5.5), near(nB, 0.2, 1),
                cbind(runif(nF, 300, 400), runif(nF, 300, 400)))
    data.frame(imageID = sprintf("img%02d", i), condition = grp, x = xy[, 1], y = xy[, 2],
               cellType = rep(c("tumour", "T", "B", "far"), c(nt, nT, nB, nF)), age = 40 + (i * 7) %% 13)
  }))
}

test_that("an adjusted design is not fitted at a radius with no variation (empty or saturated)", {
  d <- nv_cells()
  run <- function(cells, r, combine = "davies", ...) suppressMessages(spicy(cells, "condition", r = r, from = "tumour", to = "T",
    covariates = "age", adjustAbundance = TRUE, combine = combine, ...))
  for (cmb in c("davies", "maxT", "cauchy")) {
    res <- run(d, c(2, 4, 10), cmb)
    rr <- res$radiusResults
    # the unadjusted test is 0/0 at r = 2 (no T cell has a tumour neighbour) and r = 10 (every one has)
    expect_identical(rr$adjusted_for[rr$r != 4], rep("none (no_variation)", 2))
    expect_true(all(is.na(rr$p_value[rr$r != 4])))
    expect_identical(rr$adjusted_for[rr$r == 4], "abundance+covariates")
    # nothing of those radii enters the combination: as if only r = 4 were tested
    expect_identical(res$cellResults$r, 4)
    expect_identical(res$cellResults$p_value_best_radius, rr$p_value[rr$r == 4])
  }
  # the combination equals that over the radii with variation
  res <- run(d, c(2, 3, 4, 10))
  expect_identical(res$radiusResults$adjusted_for, c("none (no_variation)", rep("abundance+covariates", 2), "none (no_variation)"))
  ok <- run(d, c(3, 4))
  expect_identical(res$cellResults$p_value, ok$cellResults$p_value)
  expect_identical(res$cellResults$unadjusted_p_value, ok$cellResults$unadjusted_p_value)
  # a covariate perturbed in its last bits leaves the no-variation radii untested and the combination where it was
  for (s in 1:3) {
    set.seed(s); dp <- d; img <- unique(d$imageID)
    a <- setNames((40 + (seq_along(img) * 7) %% 13) * (1 + 1e-10 * rnorm(length(img))), img); dp$age <- unname(a[dp$imageID])
    rp <- run(dp, c(2, 3, 4, 10))
    expect_identical(rp$radiusResults$adjusted_for, res$radiusResults$adjusted_for)
    expect_equal(rp$cellResults$p_value, res$cellResults$p_value, tolerance = 1e-8)
  }
  # one radius, and the per-pair path (.cell_rows_test, as Kontextual uses it)
  for (r in c(2, 10)) {
    x <- run(d, r)$cellResults
    expect_identical(x$adjusted_for, "none (no_variation)")
    expect_true(is.na(x$p_value))
    ctx <- spicyR:::.cell_context(d, "condition", NULL, "imageID", "cellType", c("x", "y"))
    g <- spicyR:::.cell_graph(ctx, list(c("tumour", "T")), r = r, label_clustering = FALSE)
    g$image_heterogeneity <- TRUE
    Z <- matrix(d$age[ctx$first] - mean(d$age[ctx$first]), dimnames = list(NULL, "age"))
    o <- spicyR:::.cell_pair_test(ctx, g, "tumour", "T", TRUE, "cr2", TRUE, Z)
    expect_identical(o$adjusted_for, "none (no_variation)")
    expect_true(is.nan(o$test$p))
  }
})

test_that("the design test reports no variation when the response is fitted exactly", {
  set.seed(1)
  n <- 12; rows <- list(img = 0:(n - 1), O = rep(0, n), E = runif(n, 1, 3), n = runif(n, 1, 3), v = runif(n, 0.5, 1))
  rows$O <- rows$E - rows$n   # (O - E) / n = -1 in every image: the condition means fit it exactly
  g <- rep(0:1, each = n / 2); Z <- cbind(g == 0, g == 1, rnorm(n)) * 1
  a <- spicyR:::stats_design_tests(rows, 0:(n - 1), n, Z, rbind(c(-1, 1, 0)), 0, FALSE)
  expect_false(a[[1]]$ok); expect_identical(a[[1]]$reason, "no_variation")
  rows$O <- rows$O + rnorm(n, 0, 1e-6)   # a little variation is tested
  expect_true(spicyR:::stats_design_tests(rows, 0:(n - 1), n, Z, rbind(c(-1, 1, 0)), 0, FALSE)[[1]]$ok)
})

test_that("more than two conditions and survival do not test a pair with no variation", {
  d <- nv_cells(24)
  d$c3 <- c("A", "B", "C")[(as.integer(factor(d$imageID)) %% 3) + 1]
  ctx <- spicyR:::.cell_context(d, "c3", NULL, "imageID", "cellType", c("x", "y"))
  for (r in c(2, 10)) {
    g <- spicyR:::.cell_graph(ctx, list(c("tumour", "T")), r = r, label_clustering = FALSE)
    g$image_heterogeneity <- TRUE
    o <- spicyR:::.cell_pair_test(ctx, g, "tumour", "T", TRUE, "cr2", TRUE)
    expect_false(o$ok); expect_identical(o$reason, "no_variation")
    rows <- spicyR:::.cell_rows(ctx, g, "tumour", "T")
    set.seed(2); m <- length(ctx$image_labels)
    s <- spicyR:::stats_survival_test(rows, rows$unit, m, rnorm(m), rexp(m), rbinom(m, 1, 0.7), numeric(0))
    expect_false(s$ok); expect_identical(s$reason, "no_variation")
  }
  # at r = 4 the pair varies and is tested
  g <- spicyR:::.cell_graph(ctx, list(c("tumour", "T")), r = 4, label_clustering = FALSE)
  g$image_heterogeneity <- TRUE
  expect_true(spicyR:::.cell_pair_test(ctx, g, "tumour", "T", TRUE, "cr2", TRUE)$ok)
})

## tumour -> T avoids, and at r = 4 no T cell of condition A has a tumour cell within r (every A image at -1): A has no
## informative patient there. Condition B has some T cells within 4; at r = 4.8 and 5.2 some A T cells are within r.
ni_cells <- function(n_img = 16, seed = 5, levels = c("A", "B")) {
  set.seed(seed)
  do.call(rbind, lapply(seq_len(n_img), function(i) {
    g <- expand.grid(x = seq(10, 190, 20), y = seq(10, 190, 20)); nt <- sample(15:25, 1)
    g <- g[sample(nrow(g), nt), ]; tx <- g$x; ty <- g$y
    near <- function(n, lo, hi) { k <- sample(nt, n, TRUE); a <- runif(n, 0, 2 * pi); d <- runif(n, lo, hi)
      cbind(tx[k] + d * cos(a), ty[k] + d * sin(a)) }
    grp <- levels[(i - 1) %% length(levels) + 1]
    nT <- sample(30:50, 1); nc <- if (grp == "A") 0 else round(nT * runif(1, 0.05, 0.3))
    nB <- sample(40:60, 1); nF <- sample(30:50, 1)
    xy <- rbind(cbind(tx, ty), near(nc, 2.5, 3.5), near(nT - nc, 4.5, 5.5), near(nB, 0.2, 1),
                cbind(runif(nF, 300, 400), runif(nF, 300, 400)))
    data.frame(imageID = sprintf("img%02d", i), condition = grp, x = xy[, 1], y = xy[, 2],
               cellType = rep(c("tumour", "T", "B", "far"), c(nt, nT, nB, nF)), age = 40 + (i * 7) %% 13)
  }))
}

test_that("a radius where a condition has no informative patient is not tested", {
  d <- ni_cells()
  run <- function(r, ...) suppressMessages(spicy(d, "condition", r = r, from = "tumour", to = "T", ...))
  # the rows at r = 4: every image of A at -1, so A has no informative patient (B has all 8)
  ctx <- spicyR:::.cell_context(d, "condition", NULL, "imageID", "cellType", c("x", "y"))
  g <- spicyR:::.cell_graph(ctx, list(c("tumour", "T")), r = 4, label_clustering = FALSE); g$image_heterogeneity <- TRUE
  rows <- spicyR:::.cell_rows(ctx, g, "tumour", "T")
  expect_identical(attr(rows, "side"), "avoid")
  expect_true(all(rows$O[rows$group == 0L] == -1))
  expect_identical(spicyR:::stats_informative_units(rows, rows$unit, rows$group, 2L), c(0L, 8L))
  u <- spicyR:::stats_excess_test(rows, rows$unit, rows$group, length(ctx$unit_labels), TRUE, "cr2")
  expect_true(u$ok && u$se > 0 && u$coef_ref == -1)
  # one radius: not tested, unadjusted and adjusted
  x <- run(4)$cellResults
  expect_true(is.na(x$p_value)); expect_identical(x$excess_ref, -1); expect_identical(x$se, u$se)
  x <- run(4, covariates = "age", adjustAbundance = TRUE)$cellResults
  expect_identical(x$adjusted_for, "none (no_information)")
  expect_true(is.na(x$p_value)); expect_true(is.na(x$unadjusted_p_value))
  o <- spicyR:::.cell_pair_test(ctx, g, "tumour", "T", TRUE, "cr2", FALSE)
  expect_true(o$ok); expect_true(is.nan(o$test$p))
  Z <- matrix(d$age[ctx$first] - mean(d$age[ctx$first]), dimnames = list(NULL, "age"))
  o <- spicyR:::.cell_pair_test(ctx, g, "tumour", "T", TRUE, "cr2", TRUE, Z)
  expect_identical(o$adjusted_for, "none (no_information)"); expect_true(is.nan(o$test$p))
  # several radii: r = 4 is left out of the combination
  for (cmb in c("davies", "maxT", "cauchy")) for (adj in c(FALSE, TRUE)) {
    a <- run(c(4, 4.8, 5.2), combine = cmb, adjustAbundance = adj)
    b <- run(c(4.8, 5.2), combine = cmb, adjustAbundance = adj)
    rr <- a$radiusResults
    expect_true(is.na(rr$p_value[rr$r == 4])); expect_false(anyNA(rr$p_value[rr$r != 4]))
    if (adj) expect_identical(rr$adjusted_for[rr$r == 4], "none (no_information)")
    expect_identical(a$cellResults$p_value, b$cellResults$p_value)
    expect_identical(a$cellResults$r, b$cellResults$r)
    if (adj) expect_identical(a$cellResults$unadjusted_p_value, b$cellResults$unadjusted_p_value)
  }
  # the count effect has no bound and is tested as before (here O = 0 gives -E / n, which varies)
  for (r in c(4, 4.8)) {
    gc <- spicyR:::.cell_graph(ctx, list(c("tumour", "T")), r = r, label_clustering = FALSE, effect = "count")
    gc$image_heterogeneity <- TRUE
    rc <- spicyR:::.cell_rows(ctx, gc, "tumour", "T")
    uc <- spicyR:::stats_excess_test(rc, rc$unit, rc$group, length(ctx$unit_labels), TRUE, "cr2")
    x <- run(r, effect = "count")$cellResults
    expect_identical(x$p_value, uc$p); expect_false(is.na(x$p_value))
  }
  # more than two conditions: a pair where any condition has no informative patient is not tested
  d3 <- ni_cells(24, levels = c("A", "B", "C"))
  ctx3 <- spicyR:::.cell_context(d3, "condition", NULL, "imageID", "cellType", c("x", "y"))
  for (r in c(4, 4.8)) {
    g3 <- spicyR:::.cell_graph(ctx3, list(c("tumour", "T")), r = r, label_clustering = FALSE); g3$image_heterogeneity <- TRUE
    o <- spicyR:::.cell_pair_test(ctx3, g3, "tumour", "T", TRUE, "cr2", FALSE)
    if (r == 4) { expect_false(o$ok); expect_identical(o$reason, "no_information") } else expect_true(o$ok)
  }
})
