## Small synthetic checks on the R front end. Agreement with the C++ core's
## numbers is covered by the cross-check against the Python front end, which
## calls the same core; these pin the R-side enumeration and skip rules.

make_cells <- function(n_img = 8, per_type = 60, types = c("A", "B", "C"), seed = 1) {
  set.seed(seed)
  do.call(rbind, lapply(seq_len(n_img), function(i) {
    tt <- if (i == 1) types[types != "C"] else types
    n <- per_type * length(tt)
    data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000),
               cellType = rep(tt, each = per_type),
               imageID = sprintf("i%d", i),
               response = if (i %% 2 == 1) "beta" else "alpha",
               stringsAsFactors = FALSE)
  }))
}

test_that("poisson fits every unordered pair once, plus the self-pairs", {
  out <- spicy_glm(make_cells(), condition = "response", r = 60)
  key <- paste(pmin(out$results$from, out$results$to), pmax(out$results$from, out$results$to))
  expect_equal(nrow(out$results) + nrow(out$skipped), 6L)   # choose(3,2) + 3
  expect_false(any(duplicated(key)))
})

test_that("poisson is direction-invariant", {
  ab <- spicy_glm(make_cells(), condition = "response", r = 60, from = "A", to = "B")$results
  ba <- spicy_glm(make_cells(), condition = "response", r = 60, from = "B", to = "A")$results
  expect_equal(ab$log_rate_ratio, ba$log_rate_ratio, tolerance = 1e-12)
  expect_equal(ab$p_value, ba$p_value, tolerance = 1e-10)
})

test_that("binomial fits every ordered pair, and the two directions differ", {
  out <- spicy_glm(make_cells(), condition = "response", family = "binomial", k = 10)
  key <- paste(out$results$from, out$results$to, sep = "->")
  expect_setequal(key, as.vector(outer(c("A", "B", "C"), c("A", "B", "C"), paste, sep = "->")))
  expect_gt(abs(out$results$log_odds_ratio[key == "A->B"] - out$results$log_odds_ratio[key == "B->A"]), 1e-3)
})

test_that("binomial all-pairs rows equal the single-pair fits in each direction", {
  out <- spicy_glm(make_cells(), condition = "response", family = "binomial", k = 10)$results
  for (p in list(c("A", "B"), c("B", "A"))) {
    one <- spicy_glm(make_cells(), condition = "response", family = "binomial", k = 10,
                     from = p[1], to = p[2])$results
    row <- out[out$from == p[1] & out$to == p[2], ]
    expect_equal(row$log_odds_ratio, one$log_odds_ratio, tolerance = 1e-12)
    expect_equal(row$p_value, one$p_value, tolerance = 1e-12)
  }
})

test_that("binomial from and to vectors fit exactly from x to", {
  out <- spicy_glm(make_cells(), condition = "response", family = "binomial", k = 10,
                   from = "A", to = c("B", "C"))
  expect_setequal(paste(out$results$from, out$results$to, sep = "->"), c("A->B", "A->C"))
  only_to <- spicy_glm(make_cells(), condition = "response", family = "binomial", k = 10, to = "A")
  expect_setequal(paste(only_to$results$from, only_to$results$to, sep = "->"), c("A->A", "B->A", "C->A"))
})

test_that("a factor condition sets the reference from its level order", {
  cells <- make_cells()
  cells$response <- factor(cells$response, levels = c("beta", "alpha"))
  out <- spicy_glm(cells, condition = "response", r = 60)
  expect_true(all(out$results$condition_ref == "beta"))
  expect_true(all(out$results$condition_comp == "alpha"))
})

test_that("a character condition uses the sorted first level as reference", {
  out <- spicy_glm(make_cells(), condition = "response", r = 60)
  expect_true(all(out$results$condition_ref == "alpha"))
})

test_that("from and to select a single pair", {
  out <- spicy_glm(make_cells(), condition = "response", r = 60, from = "A", to = "B")
  expect_equal(nrow(out$results), 1L)
  expect_equal(out$results$from, "A")
  expect_equal(out$results$to, "B")
})

test_that("the binomial family needs k and the poisson family needs r", {
  cells <- make_cells()
  expect_error(spicy_glm(cells, condition = "response"), "positive radius")
  expect_error(spicy_glm(cells, condition = "response", family = "binomial"), "positive integer")
})

test_that("more than two conditions is an error", {
  cells <- make_cells()
  cells$response[cells$imageID == "i3"] <- "gamma"
  expect_error(spicy_glm(cells, condition = "response", r = 60), "exactly two conditions")
})

test_that("a condition that varies within an image is an error", {
  cells <- make_cells()
  cells$response[1] <- "alpha"; cells$response[2] <- "beta"
  cells$imageID[1:2] <- "i1"
  expect_error(spicy_glm(cells, condition = "response", r = 60), "constant within each image")
})

test_that("results carry the effect columns for the family and are BH ordered", {
  out <- spicy_glm(make_cells(), condition = "response", r = 60)
  expect_true(all(c("log_rate_ratio", "rate_ratio", "p_adj") %in% names(out$results)))
  expect_false(is.unsorted(out$results$p_adj))
  b <- spicy_glm(make_cells(), condition = "response", family = "binomial", k = 10)
  expect_true(all(c("log_odds_ratio", "odds_ratio") %in% names(b$results)))
})

test_that("diagnostics are returned with the documented shape", {
  out <- spicy_glm(make_cells(), condition = "response", r = 60, compute_diagnostics = TRUE)
  d <- out$diagnostics
  expect_named(d, c("pair", "patient", "image", "cross_pair"))
  expect_named(d$cross_pair, c("patient", "image"))
  expect_equal(nrow(d$pair), nrow(out$results))
  expect_true(all(c("nu", "max_influence", "patient_with_max_influence") %in% names(d$pair)))
  expect_true(all(c("l_i", "influence_i", "percentile_rank_leverage") %in% names(d$patient)))
  expect_true(all(c("l_ij", "influence_ij", "e_ij_share_within_patient") %in% names(d$image)))
  expect_true(all(c("n_pairs_present", "wilson_lower_influence") %in% names(d$cross_pair$patient)))
})

test_that("percentile ranks lie in the unit interval and leverage sums to one per group", {
  out <- spicy_glm(make_cells(), condition = "response", r = 60, compute_diagnostics = TRUE)
  pr <- out$diagnostics$patient$percentile_rank_influence
  expect_true(all(pr >= 0 & pr <= 1, na.rm = TRUE))
  by_pair_group <- split(out$diagnostics$patient$l_i,
                         paste(out$diagnostics$patient$from, out$diagnostics$patient$to,
                               out$diagnostics$patient$group))
  expect_true(all(abs(vapply(by_pair_group, sum, numeric(1)) - 1) < 1e-8))
})

test_that("diagnostics are refused outside poisson/firth/fast", {
  cells <- make_cells()
  expect_warning(spicy_glm(cells, condition = "response", family = "binomial", k = 10,
                           compute_diagnostics = TRUE), "compute_diagnostics requires")
  expect_warning(out <- spicy_glm(cells, condition = "response", r = 60, cr2_method = "naive",
                                  compute_diagnostics = TRUE), "compute_diagnostics requires")
  expect_null(out$diagnostics)
})

test_that("the Wilson interval matches the closed form", {
  w <- spicyglm:::wilson_interval(3, 20)
  z <- stats::qnorm(0.975); n <- 20; p <- 3 / 20
  lo <- (p + z^2/(2*n) - z*sqrt((p*(1-p) + z^2/(4*n))/n)) / (1 + z^2/n)
  expect_equal(w$lower, lo, tolerance = 1e-12)
})
