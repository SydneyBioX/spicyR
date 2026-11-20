## testCorrelation(): the correlation between the tests of a fit, rebuilt from pairwiseAssoc, imageWeights and
## conditionInformation, against the core's own per-patient quantities (excess_test's influence, unit_info, tau2).

cells <- sim_cells(seed = 3, n_pat = 7)
res <- spicy(cells, condition = "condition", subject = "patient", r = 30)

## Every pair's two-condition test through the internal steps, as spicy() fits it, keeping the core's per-patient
## quantities
core_fits <- function(cells, res, effect = "allocation") {
  ctx <- spicyR:::.cell_context(cells, "condition", "patient", "imageID", "cellType", c("x", "y"))
  pairs <- unname(Map(c, res$cellResults$from, res$cellResults$to))
  g <- spicyR:::.cell_graph(ctx, pairs, r = 30, label_clustering = FALSE, effect = effect)
  g$image_heterogeneity <- TRUE
  fits <- lapply(pairs, function(p) spicyR:::.cell_pair_test(ctx, g, p[1], p[2], TRUE, "cr2", FALSE))
  list(ctx = ctx, fits = fits)
}

test_that("the CR2 correlation is that of the core's influences, whose squares sum to se^2", {
  cf <- core_fits(cells, res)
  tested <- !is.na(res$cellResults$p_value)
  psi <- t(vapply(cf$fits, function(o) as.numeric(o$test$influence), numeric(length(cf$ctx$unit_labels))))
  expect_equal(sqrt(rowSums(psi^2))[tested], res$cellResults$se[tested], tolerance = 1e-12)
  u <- psi / sqrt(rowSums(psi^2))
  R0 <- tcrossprod(u[tested, ]); diag(R0) <- 1
  R <- testCorrelation(res)
  expect_equal(unname(R[tested, tested]), R0, tolerance = 1e-12)
  expect_true(all(is.na(R[!tested, ])) && all(is.na(R[, !tested])))
  expect_identical(dimnames(R), list(rownames(res$cellResults), rownames(res$cellResults)))
  expect_true(isSymmetric(R))
  expect_true(all(diag(R)[tested] == 1))
  expect_true(all(abs(R[tested, tested]) <= 1 + 1e-12))
})

test_that("conditionInformation gives each patient's weight, J / (1 + tau2 J)", {
  cf <- core_fits(cells, res)
  unit <- res$subject
  for (k in which(!is.na(res$cellResults$p_value))[1:5]) {
    o <- cf$fits[[k]]
    W <- o$test$unit_info / (1 + o$test$tau2 * o$test$unit_info)
    w <- res$imageWeights[[k]] * res$conditionInformation[k, as.character(res$condition)]
    rebuilt <- tapply(w, factor(unit, levels = cf$ctx$unit_labels), sum, na.rm = TRUE)
    expect_equal(as.numeric(rebuilt), W, tolerance = 1e-12)
  }
  expect_identical(colnames(res$conditionInformation), levels(res$condition))
  expect_identical(rownames(res$conditionInformation), res$comparisons$labels)
})

test_that("the permutation recomputes each test exactly as spicy() does", {
  tab <- res$cellResults; lab <- res$comparisons$labels; tested <- !is.na(tab$p_value)
  unit <- match(res$subject, unique(res$subject)) - 1L
  y <- do.call(cbind, unname(res$pairwiseAssoc[lab])); w <- do.call(cbind, unname(res$imageWeights[lab]))
  y[, !tested] <- NA; w[, !tested] <- NA
  obs <- as.integer(res$condition)[match(seq_along(unique(unit)) - 1L, unit)] - 1L
  x <- spicyR:::stats_test_correlation(y, w, unit, as.integer(res$condition) - 1L, length(unique(unit)),
                                       t(res$conditionInformation[lab, ]), 2L, cbind(obs, obs))
  z <- -stats::qnorm(tab$p_value / 2) * sign(tab$excess_difference)
  expect_equal(x$z[tested, 1], z[tested], tolerance = 1e-10)
  # the relabellings come from a fixed seed: the same matrix every time, with a unit diagonal
  R1 <- testCorrelation(res, "permutation", nPermutations = 100)
  expect_identical(R1, testCorrelation(res, "permutation", nPermutations = 100, cores = 2))
  expect_true(isSymmetric(R1) && all(diag(R1)[tested] == 1))
  expect_true(all(is.na(R1[!tested, ])))
})

test_that("cellWeights, the count effect, k and Kontextual give the correlation too", {
  M <- outer(cells$cellType, sort(unique(cells$cellType)), "==") * 1
  colnames(M) <- sort(unique(cells$cellType))
  M <- cbind(M, lymphocyte = M[, "T"] + M[, "B"])
  rw <- spicy(cells, condition = "condition", subject = "patient", r = 30, cellWeights = M,
              from = c("tumour", "lymphocyte"), to = c("T", "B", "macro"))
  expect_equal(dim(testCorrelation(rw)), c(nrow(rw$cellResults), nrow(rw$cellResults)))
  # with 0/1 weights a pair's test, and so its row of the correlation, is that of the cell types
  same <- intersect(rownames(rw$cellResults), rownames(res$cellResults))
  expect_equal(testCorrelation(rw)[same, same], testCorrelation(res)[same, same], tolerance = 1e-10)
  expect_silent(testCorrelation(spicy(cells, condition = "condition", subject = "patient", r = 30, effect = "count")))
  expect_silent(testCorrelation(spicy(cells, condition = "condition", subject = "patient", k = 10)))
  expect_silent(testCorrelation(spicy(cells, condition = "condition", r = 30, from = "tumour")))
  parentDf <- data.frame(from = "tumour", to = c("T", "B"), parent_name = "lymphocytes")
  parentDf$parent <- list(c("T", "B"), c("T", "B"))
  k <- kontextualEngine(cells, parentDf, condition = "condition", subject = "patient", r = 30)
  expect_equal(dim(testCorrelation(k)), c(2L, 2L))
})

test_that("testCorrelation() refuses fits whose tests it cannot rebuild", {
  expect_error(testCorrelation(spicy(cells, condition = "condition", subject = "patient", r = 30, from = "tumour",
                                     adjustAbundance = TRUE)), "unadjusted")
  expect_error(testCorrelation(spicy(cells, condition = "condition", subject = "patient", r = c(20, 30),
                                     from = "tumour")), "one radius")
  expect_error(testCorrelation(spicy(cells, condition = "condition", subject = "patient", r = 30, from = "tumour",
                                     variance = "hartung_knapp")), "cr2")
  cells$stage <- c("I", "II", "III")[match(cells$patient, unique(cells$patient)) %% 3 + 1]
  r3 <- spicy(cells, condition = "stage", subject = "patient", r = 30, from = "tumour", to = c("T", "B"))
  expect_null(r3$conditionInformation)
  expect_error(testCorrelation(r3), "two conditions")
  expect_error(testCorrelation(res, "permutation", nPermutations = 1), "nPermutations")
  expect_silent(testCorrelation(res, nPermutations = NULL))   # not used by "cr2"
  expect_error(testCorrelation(res, cores = 0), "cores")
  parentDf <- data.frame(from = "tumour", to = c("T", "B"), parent_name = "lymphocytes")
  parentDf$parent <- list(c("T", "B"), c("T", "B"))
  kh <- kontextualEngine(cells, parentDf, condition = "condition", subject = "patient", r = 30,
                         variance = "hartung_knapp")
  expect_identical(kh$variance, "hartung_knapp")
  expect_error(testCorrelation(kh), "cr2")
  # imageWeights that are not those of the tests (here: of another fit)
  bad <- res; bad$imageWeights <- spicy(cells, condition = "condition", subject = "patient", r = 20)$imageWeights
  expect_error(testCorrelation(bad), "do not give its tests")
})

test_that("patients with a missing condition take no part", {
  cm <- cells; cm$condition <- as.character(cm$condition)
  cm$condition[cm$patient %in% c("A02", "B03")] <- NA
  rm <- spicy(cm, condition = "condition", subject = "patient", r = 30, from = "tumour")
  tested <- !is.na(rm$cellResults$p_value)
  for (m in c("cr2", "permutation")) {
    R <- testCorrelation(rm, m, nPermutations = 100)
    expect_true(all(diag(R)[tested] == 1) && isSymmetric(R))
  }
})
