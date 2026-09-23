## Tests for the Binomial (fixed-k nearest-neighbour) spicyGLM backend and the
## working-variance refactor of the fast CR2 machinery.
##
## Per project convention only the default binomial path is exercised end to end
## (cr2Method = "fast", estimator = "firth" -> brglm2, fastMethod = "direct").

skip_if_not_installed("SpatialExperiment")
skip_if_not_installed("brglm2")

suppressWarnings(suppressPackageStartupMessages(library(SpatialExperiment)))

## --- synthetic data ---------------------------------------------------------

make_spe <- function(seed = 42, n_per_cond = 6, cells_per_img = 400,
                     enrich_cond = "B") {
  set.seed(seed)
  conds <- c("A", "B")
  blocks <- list()
  for (cond in conds) {
    for (s in seq_len(n_per_cond)) {
      sid <- paste0(cond, "_s", s)
      nc <- cells_per_img
      ct <- sample(c("Tcell", "Tumour", "Bcell"), nc, replace = TRUE,
                   prob = c(0.3, 0.4, 0.3))
      x <- runif(nc, 0, 1000); y <- runif(nc, 0, 1000)
      if (cond == enrich_cond) {
        tum <- ct == "Tumour"
        cx <- mean(x[tum]); cy <- mean(y[tum])
        tt <- which(ct == "Tcell")
        mv <- sample(tt, length(tt) %/% 2)
        x[mv] <- 0.5 * x[mv] + 0.5 * cx
        y[mv] <- 0.5 * y[mv] + 0.5 * cy
      }
      blocks[[length(blocks) + 1]] <- data.frame(
        imageID = sid, subject = sid, condition = cond,
        cellType = ct, x = x, y = y, stringsAsFactors = FALSE
      )
    }
  }
  df <- do.call(rbind, blocks)
  SpatialExperiment(
    assays = list(counts = matrix(0, 1, nrow(df))),
    colData = DataFrame(df[, c("imageID", "subject", "condition", "cellType")]),
    spatialCoords = as.matrix(df[, c("x", "y")])
  )
}

run_glm <- function(spe, ..., quiet = TRUE) {
  f <- function() spicyGLM(spe, condition = "condition", subject = "subject",
                           imageID = "imageID", cellType = "cellType",
                           spatialCoords = c("x", "y"), cores = 1, ...)
  if (quiet) suppressMessages(suppressWarnings(f())) else f()
}

## --- Poisson: refactor is behaviour-preserving -----------------------------

test_that("Poisson path: direct and dpr1 give identical CR2 results", {
  spe <- make_spe()
  d <- run_glm(spe, r = 50, from = "Tcell", to = "Tumour",
               family = "poisson", fastMethod = "direct")
  p <- run_glm(spe, r = 50, from = "Tcell", to = "Tumour",
               family = "poisson", fastMethod = "dpr1")
  expect_equal(d$GLMresults$logRateRatio, p$GLMresults$logRateRatio, tolerance = 1e-10)
  expect_equal(d$GLMresults$p.value,      p$GLMresults$p.value,      tolerance = 1e-10)
})

test_that("Poisson result frame keeps its schema and gains a family column", {
  spe <- make_spe()
  res <- run_glm(spe, r = 50, from = "Tcell", to = "Tumour", family = "poisson")
  expect_identical(res$family, "poisson")
  expect_true(all(c("logRateRatio", "rateRatio", "family") %in% names(res$GLMresults)))
  expect_false(any(c("logOddsRatio", "oddsRatio") %in% names(res$GLMresults)))
  expect_identical(unique(res$GLMresults$family), "poisson")
})

## --- Binomial: happy path -------------------------------------------------

test_that("Binomial single pair fits and reports an odds ratio", {
  spe <- make_spe()
  res <- run_glm(spe, family = "binomial", k = 10, from = "Tcell", to = "Tumour")
  expect_identical(res$family, "binomial")
  expect_true(all(c("logOddsRatio", "oddsRatio", "family") %in% names(res$GLMresults)))
  expect_false(any(c("logRateRatio", "rateRatio") %in% names(res$GLMresults)))
  expect_equal(res$GLMresults$oddsRatio, exp(res$GLMresults$logOddsRatio))
  expect_gt(res$GLMresults$p.value, 0)
  expect_lt(res$GLMresults$p.value, 1)
  expect_identical(res$GLMresults$estimator, "firth")
})

test_that("Binomial all-pairs run produces one row per ordered pair", {
  spe <- make_spe()
  res <- run_glm(spe, family = "binomial", k = 12)
  expect_equal(nrow(res$GLMresults) + NROW(res$skipped), 9L)
  expect_true(all(res$GLMresults$family == "binomial"))
  expect_true("p.adj" %in% names(res$GLMresults))
  expect_null(res$diagnostics)
})

## --- directionality ----------------------------------------------------------

test_that("Binomial fits both directions and they differ", {
  spe <- make_spe()
  res <- run_glm(spe, family = "binomial", k = 12)$GLMresults
  key <- paste(res$from, res$to, sep = "->")
  expect_setequal(key, as.vector(outer(c("Tcell", "Tumour", "Bcell"),
                                       c("Tcell", "Tumour", "Bcell"), paste, sep = "->")))
  ab <- res$logOddsRatio[key == "Tcell->Tumour"]
  ba <- res$logOddsRatio[key == "Tumour->Tcell"]
  expect_gt(abs(ab - ba), 1e-3)
})

test_that("Binomial all-pairs rows equal the single-pair fits in each direction", {
  spe <- make_spe()
  res <- run_glm(spe, family = "binomial", k = 12)$GLMresults
  for (p in list(c("Tcell", "Tumour"), c("Tumour", "Tcell"))) {
    one <- run_glm(spe, family = "binomial", k = 12, from = p[1], to = p[2])$GLMresults
    row <- res[res$from == p[1] & res$to == p[2], ]
    expect_equal(row$logOddsRatio, one$logOddsRatio, tolerance = 1e-10)
    expect_equal(row$p.value, one$p.value, tolerance = 1e-10)
  }
})

test_that("Binomial from/to vectors fit exactly expand.grid(from, to)", {
  spe <- make_spe()
  res <- run_glm(spe, family = "binomial", k = 12, from = "Tcell",
                 to = c("Tumour", "Bcell"))$GLMresults
  expect_setequal(paste(res$from, res$to, sep = "->"), c("Tcell->Tumour", "Tcell->Bcell"))
})

test_that("Poisson keeps one row per unordered pair and is direction-invariant", {
  spe <- make_spe()
  res <- run_glm(spe, family = "poisson", r = 50)
  expect_equal(nrow(res$GLMresults) + NROW(res$skipped), 6L)
  ab <- run_glm(spe, family = "poisson", r = 50, from = "Tcell", to = "Tumour")$GLMresults
  ba <- run_glm(spe, family = "poisson", r = 50, from = "Tumour", to = "Tcell")$GLMresults
  expect_equal(ab$logRateRatio, ba$logRateRatio, tolerance = 1e-10)
  expect_equal(ab$p.value, ba$p.value, tolerance = 1e-8)
})

test_that("Binomial direct and dpr1 agree", {
  spe <- make_spe()
  a <- run_glm(spe, family = "binomial", k = 10, from = "Tcell", to = "Tumour",
               fastMethod = "direct")
  b <- run_glm(spe, family = "binomial", k = 10, from = "Tcell", to = "Tumour",
               fastMethod = "dpr1")
  expect_equal(a$GLMresults$logOddsRatio, b$GLMresults$logOddsRatio, tolerance = 1e-8)
  expect_equal(a$GLMresults$p.value,      b$GLMresults$p.value,      tolerance = 1e-8)
})

## --- Binomial: argument handling ----------------------------------------

test_that("family = 'binomial' requires k", {
  spe <- make_spe()
  expect_error(
    suppressMessages(spicyGLM(spe, condition = "condition", subject = "subject",
                              imageID = "imageID", cellType = "cellType",
                              spatialCoords = c("x", "y"), family = "binomial",
                              from = "Tcell", to = "Tumour", cores = 1)),
    "requires `k`"
  )
})

test_that("closed-form Firth is coerced to brglm2 for binomial", {
  spe <- make_spe()
  expect_message(
    suppressWarnings(spicyGLM(spe, condition = "condition", subject = "subject",
                              imageID = "imageID", cellType = "cellType",
                              spatialCoords = c("x", "y"), family = "binomial",
                              k = 10, from = "Tcell", to = "Tumour", cores = 1,
                              firthBackend = "closed_form")),
    "no closed form"
  )
})

test_that("computeDiagnostics is refused for binomial", {
  spe <- make_spe()
  expect_message(
    suppressWarnings(spicyGLM(spe, condition = "condition", subject = "subject",
                              imageID = "imageID", cellType = "cellType",
                              spatialCoords = c("x", "y"), family = "binomial",
                              k = 10, from = "Tcell", to = "Tumour", cores = 1,
                              computeDiagnostics = TRUE)),
    "not supported for family = 'binomial'"
  )
  res <- run_glm(spe, family = "binomial", k = 10, from = "Tcell", to = "Tumour",
                 computeDiagnostics = TRUE)
  expect_null(res$diagnostics)
})

## --- Binomial: image-level drop guards ---------------------------------

test_that("images with <= k cells are dropped with a warning", {
  spe <- make_spe()
  tiny <- data.frame(imageID = "A_tiny", subject = "A_tiny", condition = "A",
                     cellType = c("Tcell", "Tumour", "Tcell", "Bcell", "Tumour"),
                     x = runif(5, 0, 50), y = runif(5, 0, 50))
  extra <- SpatialExperiment(
    assays = list(counts = matrix(0, 1, 5)),
    colData = DataFrame(tiny[, c("imageID", "subject", "condition", "cellType")]),
    spatialCoords = as.matrix(tiny[, c("x", "y")])
  )
  spe2 <- cbind(spe, extra)
  expect_warning(
    suppressMessages(spicyGLM(spe2, condition = "condition", subject = "subject",
                              imageID = "imageID", cellType = "cellType",
                              spatialCoords = c("x", "y"), family = "binomial",
                              k = 10, from = "Tcell", to = "Tumour", cores = 1)),
    "A_tiny"
  )
})

test_that("a single-cell-type image is dropped for a self-pair (p0 = 1)", {
  spe <- make_spe()
  mono <- data.frame(imageID = "B_mono", subject = "B_mono", condition = "B",
                     cellType = rep("Tumour", 60),
                     x = runif(60, 0, 100), y = runif(60, 0, 100))
  extra <- SpatialExperiment(
    assays = list(counts = matrix(0, 1, 60)),
    colData = DataFrame(mono[, c("imageID", "subject", "condition", "cellType")]),
    spatialCoords = as.matrix(mono[, c("x", "y")])
  )
  spe2 <- cbind(spe, extra)
  expect_warning(
    suppressMessages(spicyGLM(spe2, condition = "condition", subject = "subject",
                              imageID = "imageID", cellType = "cellType",
                              spatialCoords = c("x", "y"), family = "binomial",
                              k = 10, from = "Tumour", to = "Tumour", cores = 1)),
    "single cell type|non-finite logit"
  )
})

## --- Binomial: separation boundary detection --------------------------

bin_ctp <- data.frame(
  cellType = c("T", "Tum"), condition = rep(c("A", "B"), each = 2),
  n_images_with_type = 6, n_images_total = 6, stringsAsFactors = FALSE
)

bin_pair_df <- function(nA, nB, kk = 8, ncell = 25) {
  one <- function(sid, cond, nvec) data.frame(
    from = "T", to = "Tum", condition = cond, imageID = sid, subject = sid,
    n = nvec, k = kk, p0 = 0.3, stringsAsFactors = FALSE
  )
  rows <- c(
    lapply(1:6, function(s) one(paste0("A", s), "A", nA(ncell))),
    lapply(1:6, function(s) one(paste0("B", s), "B", nB(ncell)))
  )
  d <- do.call(rbind, rows)
  d$condition <- factor(d$condition)
  d
}

fit_bin <- function(d, estimator) {
  buildGLM(d, oneToOne = TRUE, cr2Method = "fast", fastMethod = "direct",
           estimator = estimator, cellTypePresence = bin_ctp,
           family = "binomial", k = 8)
}

test_that("floor separation in one condition -> one_condition_zero", {
  set.seed(1)
  d <- bin_pair_df(nA = function(n) rep(0L, n), nB = function(n) rpois(n, 2))
  firth <- suppressWarnings(fit_bin(d, "firth"))
  expect_identical(firth$mle_skip_reason, "one_condition_zero")
  expect_true(firth$mle_would_skip)
  expect_true(is.finite(firth$logOddsRatio))
  mle <- suppressWarnings(fit_bin(d, "mle"))
  expect_true("reason" %in% names(mle))
  expect_identical(mle$reason, "one_condition_zero")
})

test_that("ceiling separation in one condition -> one_condition_max", {
  set.seed(2)
  d <- bin_pair_df(nA = function(n) rep(8L, n), nB = function(n) rbinom(n, 8, 0.4))
  firth <- suppressWarnings(fit_bin(d, "firth"))
  expect_identical(firth$mle_skip_reason, "one_condition_max")
  expect_true(is.finite(firth$logOddsRatio))
  mle <- suppressWarnings(fit_bin(d, "mle"))
  expect_identical(mle$reason, "one_condition_max")
})

test_that("both conditions saturated -> all_max, Firth log OR ~ 0", {
  d <- bin_pair_df(nA = function(n) rep(8L, n), nB = function(n) rep(8L, n))
  firth <- suppressWarnings(fit_bin(d, "firth"))
  expect_identical(firth$mle_skip_reason, "all_max")
  expect_lt(abs(firth$logOddsRatio), 1e-6)
})

test_that("opposite boundaries -> all_boundary", {
  d <- bin_pair_df(nA = function(n) rep(0L, n), nB = function(n) rep(8L, n))
  firth <- suppressWarnings(fit_bin(d, "firth"))
  expect_identical(firth$mle_skip_reason, "all_boundary")
})

test_that("Poisson still only flags the floor (no ceiling concept)", {
  pois_ctp <- bin_ctp
  one <- function(sid, cond, nvec) data.frame(
    from = "T", to = "Tum", condition = cond, imageID = sid, subject = sid,
    n = nvec, density = 1.5, stringsAsFactors = FALSE
  )
  set.seed(5)
  rows <- c(
    lapply(1:6, function(s) one(paste0("A", s), "A", rep(0L, 25))),
    lapply(1:6, function(s) one(paste0("B", s), "B", rpois(25, 3)))
  )
  d <- do.call(rbind, rows); d$condition <- factor(d$condition)
  firth <- suppressWarnings(buildGLM(d, oneToOne = TRUE, cr2Method = "fast",
                                     fastMethod = "direct", estimator = "firth",
                                     cellTypePresence = pois_ctp, family = "poisson"))
  expect_identical(firth$mle_skip_reason, "one_condition_zero")
})

## --- fast CR2 vs clubSandwich (independent oracle for variance AND d.o.f.) --

test_that("Binomial fast CR2 matches clubSandwich on the SE and Satterthwaite df", {
  skip_if_not_installed("clubSandwich")

  ## one cell-type pair, ~14 clusters, plain binomial GLM (MLE, no Firth) so the
  ## comparison against clubSandwich is exact
  set.seed(2024)
  kk <- 10L
  rows <- lapply(seq_len(14), function(i) {
    cond <- if (i <= 7) "A" else "B"
    ncell <- sample(20:40, 1)
    p0 <- runif(1, 0.20, 0.50)
    pit <- plogis(qlogis(p0) + if (cond == "B") 0.35 else 0)
    data.frame(condition = cond, imageID = paste0("img", i),
               n = rbinom(ncell, kk, pit), k = kk, p0 = p0,
               stringsAsFactors = FALSE)
  })
  d <- do.call(rbind, rows)
  d$condition <- factor(d$condition)

  fit <- glm(cbind(n, k - n) ~ 0 + condition, offset = qlogis(p0),
             family = binomial(), data = d)
  logRR <- unname(coef(fit)[2] - coef(fit)[1])

  V_cs <- clubSandwich::vcovCR(fit, cluster = d$imageID, type = "CR2")
  se_cs <- sqrt(as.numeric(matrix(c(-1, 1), 1) %*% V_cs %*% c(-1, 1)))
  wt_cs <- clubSandwich::Wald_test(fit, matrix(c(-1, 1), 1), V_cs, tidy = TRUE)

  pats <- lapply(unique(d$imageID), build_patient, cluster_vec = d$imageID,
                 df_result = d, fit = fit, condition_levels = levels(d$condition),
                 family = "binomial", k = kk)
  V_fast <- vcovCR2_fast_multi(pats, method = "direct")
  wr <- waldTest_CR2_fast(logRR, V_fast, pats, method = "direct")

  expect_equal(sqrt(wr$v_hat), se_cs, tolerance = 1e-5)
  expect_equal(wr$df, wt_cs$df_denom, tolerance = 1e-4)   # the d.o.f. -- regression guard
  expect_equal(wr$p.value, wt_cs$p_val, tolerance = 1e-4)

  ## dpr1 route must agree with direct
  wr2 <- waldTest_CR2_fast(logRR, vcovCR2_fast_multi(pats, method = "dpr1"),
                           pats, method = "dpr1")
  expect_equal(wr2$df, wr$df, tolerance = 1e-6)
})
