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

test_that("every unordered pair is fitted once, plus the self-pairs", {
  out <- spicy_glm(make_cells(), condition = "response", r = 60)
  key <- paste(pmin(out$results$from, out$results$to), pmax(out$results$from, out$results$to))
  expect_equal(nrow(out$results) + nrow(out$skipped), 6L)   # choose(3,2) + 3
  expect_false(any(duplicated(key)))
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
