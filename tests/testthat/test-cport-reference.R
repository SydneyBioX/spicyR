## Reference fixtures + known outputs for validating a future C port of
## spicyGLM's fast CR2 machinery. Two data sources:
##  - a single deterministic draw from the production simulation DGP
##    (~/Documents/spicyGLM/sim_common.R's simulate_once()/DESIGN, "signal"
##    scenario, seed = 1, as of 2026-09) -- a synthetic dataset with a known
##    embedded effect (Group2's B-around-A smoothing sigma is inflated by
##    lambda/3 relative to Group1).
##  - diabetesData (spicyR's own bundled real dataset), subset to two `stage`
##    levels and the Tc -> Th pair (the same pair used by the existing
##    golden-snapshot test in test-spicyR.R).
##
## Every fit uses cr2Method = "fast", estimator = "firth", fastMethod =
## "direct" -- the only tested/supported production path (Poisson uses
## firthBackend = "closed_form"; Binomial is forced to "brglm2", since Firth
## has no closed form there).
##
## The hardcoded numbers below were verified once during development against
## clubSandwich::vcovCR(type="CR2")/Wald_test() on an equivalent plain MLE
## glm() fit -- vcovCR2_fast_multi()/waldTest_CR2_fast() only consume
## fitted()/residuals(), so validating them on the MLE fit carries over to the
## Firth fit used for the production literals. That oracle comparison is kept
## as a live test below; every other check is a plain
## expect_equal(x, <literal>, tolerance = ...) with no clubSandwich dependency.

skip_if_not_installed("SpatialExperiment")
skip_if_not_installed("brglm2")
skip_if_not_installed("spatstat.geom")
skip_if_not_installed("spatstat.random")
skip_if_not_installed("spatstat.explore")

suppressWarnings(suppressPackageStartupMessages(library(SpatialExperiment)))

## --- shared helper: MLE (oracle-compatible) + Firth (production) CR2 fit ---

fit_cr2_pair <- function(cells, condition, subject, imageID, cellType,
                         spatialCoords, from, to, r, window, family, k, oneToOne) {
  ctp <- computeCellTypePresence(cells, condition, imageID, cellType)
  dfPair <- suppressWarnings(modelDataGen(
    cells = cells, condition = condition, subject = subject, from = from, to = to,
    r = r, imageID = imageID, cellType = cellType, spatialCoords = spatialCoords,
    window = window, cores = 1, oneToOne = oneToOne, cellTypePresence = ctp,
    family = family, k = k
  ))
  clusterVec <- if (oneToOne) dfPair$imageID else dfPair$subject
  condLevels <- levels(droplevels(factor(dfPair$condition)))

  fit_and_wald <- function(estimator, backend) {
    fit <- fit_pair(dfPair, estimator = estimator, backend = backend, family = family)$fit
    patients <- lapply(unique(clusterVec), build_patient, cluster_vec = clusterVec,
                       df_result = dfPair, fit = fit, condition_levels = condLevels,
                       family = family, k = k)
    V <- vcovCR2_fast_multi(patients, method = "direct")
    beta <- stats::coef(fit)
    logEff <- unname(beta[2] - beta[1])
    wr <- waldTest_CR2_fast(logEff, V, patients, method = "direct")
    list(fit = fit, patients = patients, logEff = logEff,
        se = sqrt(wr$v_hat), df = wr$df, p.value = wr$p.value)
  }

  firthBackend <- if (family == "binomial") "brglm2" else "closed_form"
  list(dfPair = dfPair, clusterVec = clusterVec,
      mle = fit_and_wald("mle", "glm"), firth = fit_and_wald("firth", firthBackend))
}

cr2_eigenvalues_ok <- function(patients) {
  Sg <- compute_group_totals(patients)
  vapply(patients, function(p) {
    reduced <- reduced_dpr1(p)
    eig <- eig_Gbar_direct(reduced$Lambda_bar, reduced$f_bar, Sg[p$group])
    Abar <- assemble_Abar(p$var_hat, eig)
    n_i <- length(p$var_hat)
    all(eigen(Abar + diag(n_i))$values >= 1 - 1e-8)
  }, logical(1))
}

## --- fixture 1: simulated data (known embedded signal) ---------------------
## Copied from the production nSim=1000 simulation harness's DGP: Poisson
## cell counts per image, B cells density-modulated around A cells.

sim_simulate_once <- function(i, lambda, nPatients, nIm, counts, window,
                              scenario = c("signal", "null_equal_coloc")) {
  scenario <- match.arg(scenario)
  set.seed(i)
  g1 <- rpois(nPatients / 2, lambda)
  g2 <- if (scenario == "signal") rpois(nPatients / 2, lambda + lambda / 3) else rpois(nPatients / 2, lambda)
  adjustSigma <- c(g1, g2) + 1
  x <- y <- cellType <- imageID <- condition <- subject <- NULL
  for (p in seq_len(nPatients)) {
    for (j in seq_len(nIm)) {
      repeat { sA <- rpois(1, sample(counts, 1)); if (sA > 0) break }
      repeat { sB <- rpois(1, sample(counts, 1)); if (sB > 0) break }
      A <- spatstat.random::rpoispp(sA / spatstat.geom::area.owin(window), win = window)
      aDens <- spatstat.explore::density.ppp(A, sigma = adjustSigma[p], kernel = "disc")
      T <- spatstat.geom::integral.im(aDens)
      aDens$v <- pmax(aDens$v, 0) * (sB / T)
      B <- spatstat.random::rpoispp(aDens)
      lbl <- paste(p, j, sep = "_")
      x <- c(x, A$x, B$x); y <- c(y, A$y, B$y)
      cellType <- c(cellType, rep("A", A$n), rep("B", B$n))
      imageID  <- c(imageID, rep(lbl, A$n + B$n))
      group    <- if (p <= nPatients / 2) "Group1" else "Group2"
      condition <- c(condition, rep(group, A$n + B$n))
      subject   <- c(subject, rep(p, A$n + B$n))
    }
  }
  data.frame(x, y, cellType = factor(cellType), imageID = factor(imageID),
            condition = factor(condition), subject = factor(subject))
}

sim_build_spe <- function(cellExp) {
  counts <- matrix(0, nrow = 1, ncol = nrow(cellExp), dimnames = list("dummy_feature", NULL))
  SpatialExperiment(
    assays = list(counts = counts),
    spatialCoords = cbind(x = cellExp$x, y = cellExp$y),
    colData = cellExp[, c("cellType", "imageID", "subject", "condition")]
  )
}

SIM_DESIGN <- list(
  nPatients = 40, lambda = 40, r = 50, k = 5,
  counts = seq(20, 400, by = 10),
  window = spatstat.geom::owin(c(0, 1000), c(0, 1000))
)

sim_cellExp <- sim_simulate_once(1, lambda = SIM_DESIGN$lambda, nPatients = SIM_DESIGN$nPatients,
                                 nIm = 1, counts = SIM_DESIGN$counts, window = SIM_DESIGN$window,
                                 scenario = "signal")
sim_spe <- sim_build_spe(sim_cellExp)

run_glm_sim <- function(..., quiet = TRUE) {
  f <- function() spicyGLM(sim_spe, condition = "condition", subject = "subject",
                           imageID = "imageID", cores = 1, ...)
  if (quiet) suppressMessages(suppressWarnings(f())) else f()
}

## --- fixture 2: diabetesData subset (stage: Onset vs Long-duration) --------

data("diabetesData", package = "spicyR")
diabetes_subset <- local({
  sub <- diabetesData[, diabetesData$stage %in% c("Onset", "Long-duration")]
  sub$stage <- droplevels(factor(sub$stage))
  sub
})

run_glm_diabetes <- function(..., quiet = TRUE) {
  f <- function() spicyGLM(diabetes_subset, condition = "stage", subject = "case",
                           imageID = "imageID", cellType = "cellType",
                           spatialCoords = c("x", "y"), from = "Tc", to = "Th",
                           cores = 1, ...)
  if (quiet) suppressMessages(suppressWarnings(f())) else f()
}

## --- reproducibility guard on the simulated fixture ------------------------

test_that("the simulated fixture's DGP is deterministic given its seed", {
  d1 <- sim_simulate_once(1, lambda = SIM_DESIGN$lambda, nPatients = SIM_DESIGN$nPatients,
                          nIm = 1, counts = SIM_DESIGN$counts, window = SIM_DESIGN$window,
                          scenario = "signal")
  d2 <- sim_simulate_once(1, lambda = SIM_DESIGN$lambda, nPatients = SIM_DESIGN$nPatients,
                          nIm = 1, counts = SIM_DESIGN$counts, window = SIM_DESIGN$window,
                          scenario = "signal")
  expect_identical(d1, d2)
  expect_equal(nrow(d1), 15866)
})

## --- Poisson: simulated data ------------------------------------------------

test_that("spicyGLM Poisson on the simulated signal fixture matches known output", {
  res <- run_glm_sim(from = "B", to = "A", r = SIM_DESIGN$r, family = "poisson",
                     cr2Method = "fast", estimator = "firth", fastMethod = "direct")
  expect_equal(res$GLMresults$logRateRatio, -0.0886563689, tolerance = 1e-8)
  expect_equal(res$GLMresults$rateRatio,     0.9151600,    tolerance = 1e-6)
  expect_equal(res$GLMresults$p.value,       0.0738397252, tolerance = 1e-8)
})

test_that("fast CR2 SE/df for the simulated Poisson fit match clubSandwich (oracle)", {
  skip_if_not_installed("clubSandwich")

  fit <- fit_cr2_pair(sim_spe, condition = "condition", subject = "subject",
                      imageID = "imageID", cellType = "cellType",
                      spatialCoords = c("x", "y"), from = "B", to = "A",
                      r = SIM_DESIGN$r, window = "convex", family = "poisson",
                      k = NULL, oneToOne = TRUE)

  L <- matrix(c(-1, 1), nrow = 1)
  V_cs <- clubSandwich::vcovCR(fit$mle$fit, cluster = fit$clusterVec, type = "CR2")
  se_cs <- sqrt(as.numeric(L %*% V_cs %*% t(L)))
  wt_cs <- clubSandwich::Wald_test(fit$mle$fit, L, V_cs, tidy = TRUE)

  expect_equal(fit$mle$se, se_cs,           tolerance = 1e-5)
  expect_equal(fit$mle$df, wt_cs$df_denom,  tolerance = 1e-4)
  expect_equal(fit$mle$p.value, wt_cs$p_val, tolerance = 1e-4)

  ## having verified the machinery above (on the MLE fit), these are the
  ## production (Firth) fit's known SE/df -- what a C port should reproduce:
  expect_equal(fit$firth$se, 0.0473157016, tolerance = 1e-8)
  expect_equal(fit$firth$df, 22.805730,    tolerance = 1e-4)
})

## --- Binomial: simulated data -----------------------------------------------

test_that("spicyGLM Binomial on the simulated signal fixture matches known output", {
  res <- run_glm_sim(from = "B", to = "A", family = "binomial", k = SIM_DESIGN$k,
                     cr2Method = "fast", estimator = "firth", fastMethod = "direct")
  expect_equal(res$GLMresults$logOddsRatio, -0.0067005582, tolerance = 1e-8)
  expect_equal(res$GLMresults$oddsRatio,     0.9933218,    tolerance = 1e-6)
  expect_equal(res$GLMresults$p.value,       0.7745445391, tolerance = 1e-8)
})

test_that("fast CR2 SE/df for the simulated Binomial fit match clubSandwich (oracle)", {
  skip_if_not_installed("clubSandwich")

  fit <- fit_cr2_pair(sim_spe, condition = "condition", subject = "subject",
                      imageID = "imageID", cellType = "cellType",
                      spatialCoords = c("x", "y"), from = "B", to = "A",
                      r = NULL, window = "convex", family = "binomial",
                      k = SIM_DESIGN$k, oneToOne = TRUE)

  L <- matrix(c(-1, 1), nrow = 1)
  V_cs <- clubSandwich::vcovCR(fit$mle$fit, cluster = fit$clusterVec, type = "CR2")
  se_cs <- sqrt(as.numeric(L %*% V_cs %*% t(L)))
  wt_cs <- clubSandwich::Wald_test(fit$mle$fit, L, V_cs, tidy = TRUE)

  expect_equal(fit$mle$se, se_cs,           tolerance = 1e-5)
  expect_equal(fit$mle$df, wt_cs$df_denom,  tolerance = 1e-4)
  expect_equal(fit$mle$p.value, wt_cs$p_val, tolerance = 1e-4)

  expect_equal(fit$firth$se, 0.0231532819, tolerance = 1e-8)
  expect_equal(fit$firth$df, 26.299852,    tolerance = 1e-4)
})

## --- diagnostic identities (Poisson only; simulated fixture) ---------------

test_that("Firth closed form beta = log((Y_g+0.5)/D_g) holds on real modelData", {
  res <- run_glm_sim(from = "B", to = "A", r = SIM_DESIGN$r, family = "poisson",
                     cr2Method = "fast", estimator = "firth", fastMethod = "direct",
                     storeModelData = TRUE)
  md <- res$modelData
  Yg <- tapply(md$n, md$condition, sum)
  Dg <- tapply(md$density, md$condition, sum)
  recomputedLogRR <- unname(diff(log((Yg + 0.5) / Dg)))
  expect_equal(recomputedLogRR, res$GLMresults$logRateRatio, tolerance = 1e-10)
})

test_that("leverage l_i matches known values for representative patients", {
  res <- run_glm_sim(from = "B", to = "A", r = SIM_DESIGN$r, family = "poisson",
                     cr2Method = "fast", estimator = "firth", fastMethod = "direct",
                     computeDiagnostics = TRUE)
  patTbl <- res$diagnostics$patient
  minRow <- patTbl[which.min(patTbl$l_i), ]
  maxRow <- patTbl[which.max(patTbl$l_i), ]
  expect_identical(minRow$patient_id, "29_1")
  expect_equal(minRow$l_i, 0.001935254, tolerance = 1e-6)
  expect_identical(maxRow$patient_id, "34_1")
  expect_equal(maxRow$l_i, 0.1491479, tolerance = 1e-6)
})

test_that("CR2 A_i eigenvalues are >= 1 (never dampens residuals)", {
  fit <- fit_cr2_pair(sim_spe, condition = "condition", subject = "subject",
                      imageID = "imageID", cellType = "cellType",
                      spatialCoords = c("x", "y"), from = "B", to = "A",
                      r = SIM_DESIGN$r, window = "convex", family = "poisson",
                      k = NULL, oneToOne = TRUE)
  expect_true(all(cr2_eigenvalues_ok(fit$firth$patients)))
})

## --- Poisson: diabetesData (Tc -> Th) ---------------------------------------

test_that("spicyGLM Poisson on diabetesData (Tc->Th) matches known output", {
  res <- run_glm_diabetes(r = 50, family = "poisson",
                          cr2Method = "fast", estimator = "firth", fastMethod = "direct")
  expect_identical(res$GLMresults$conditionRef, "Onset")
  expect_identical(res$GLMresults$conditionComp, "Long-duration")
  expect_equal(res$GLMresults$logRateRatio, 0.3071080276, tolerance = 1e-8)
  expect_equal(res$GLMresults$rateRatio,    1.359488,     tolerance = 1e-6)
  expect_equal(res$GLMresults$p.value,      0.2832142993, tolerance = 1e-8)
})

test_that("fast CR2 SE/df for diabetesData Poisson fit match clubSandwich (oracle)", {
  skip_if_not_installed("clubSandwich")

  fit <- fit_cr2_pair(diabetes_subset, condition = "stage", subject = "case",
                      imageID = "imageID", cellType = "cellType",
                      spatialCoords = c("x", "y"), from = "Tc", to = "Th",
                      r = 50, window = "convex", family = "poisson",
                      k = NULL, oneToOne = FALSE)

  L <- matrix(c(-1, 1), nrow = 1)
  V_cs <- clubSandwich::vcovCR(fit$mle$fit, cluster = fit$clusterVec, type = "CR2")
  se_cs <- sqrt(as.numeric(L %*% V_cs %*% t(L)))
  wt_cs <- clubSandwich::Wald_test(fit$mle$fit, L, V_cs, tidy = TRUE)

  ## only 8 clusters (4 per group) here, vs. 40 in the simulated fixture, so
  ## the fast/clubSandwich agreement is a little looser numerically
  expect_equal(fit$mle$se, se_cs,           tolerance = 1e-4)
  expect_equal(fit$mle$df, wt_cs$df_denom,  tolerance = 1e-3)
  expect_equal(fit$mle$p.value, wt_cs$p_val, tolerance = 1e-3)

  expect_equal(fit$firth$se, 0.2359375684, tolerance = 1e-6)
  expect_equal(fit$firth$df, 3.027978,     tolerance = 1e-3)
})

## --- Binomial: diabetesData (Tc -> Th) --------------------------------------

test_that("spicyGLM Binomial on diabetesData (Tc->Th) matches known output", {
  res <- run_glm_diabetes(family = "binomial", k = 15,
                          cr2Method = "fast", estimator = "firth", fastMethod = "direct")
  expect_equal(res$GLMresults$logOddsRatio, 0.4066833365, tolerance = 1e-8)
  expect_equal(res$GLMresults$oddsRatio,    1.501828,     tolerance = 1e-6)
  expect_equal(res$GLMresults$p.value,      0.1457285988, tolerance = 1e-8)
})

test_that("fast CR2 SE/df for diabetesData Binomial fit match clubSandwich (oracle)", {
  skip_if_not_installed("clubSandwich")

  fit <- fit_cr2_pair(diabetes_subset, condition = "stage", subject = "case",
                      imageID = "imageID", cellType = "cellType",
                      spatialCoords = c("x", "y"), from = "Tc", to = "Th",
                      r = NULL, window = "convex", family = "binomial",
                      k = 15, oneToOne = FALSE)

  L <- matrix(c(-1, 1), nrow = 1)
  V_cs <- clubSandwich::vcovCR(fit$mle$fit, cluster = fit$clusterVec, type = "CR2")
  se_cs <- sqrt(as.numeric(L %*% V_cs %*% t(L)))
  wt_cs <- clubSandwich::Wald_test(fit$mle$fit, L, V_cs, tidy = TRUE)

  expect_equal(fit$mle$se, se_cs,           tolerance = 1e-4)
  expect_equal(fit$mle$df, wt_cs$df_denom,  tolerance = 1e-3)
  expect_equal(fit$mle$p.value, wt_cs$p_val, tolerance = 1e-3)

  expect_equal(fit$firth$se, 0.2155877761, tolerance = 1e-6)
  expect_equal(fit$firth$df, 3.361842,     tolerance = 1e-3)
})

test_that("Firth closed form holds on diabetesData modelData", {
  res <- run_glm_diabetes(r = 50, family = "poisson",
                          cr2Method = "fast", estimator = "firth", fastMethod = "direct",
                          storeModelData = TRUE)
  md <- res$modelData
  Yg <- tapply(md$n, md$condition, sum)
  Dg <- tapply(md$density, md$condition, sum)
  recomputedLogRR <- unname(diff(log((Yg + 0.5) / Dg)))
  expect_equal(recomputedLogRR, res$GLMresults$logRateRatio, tolerance = 1e-10)
})
