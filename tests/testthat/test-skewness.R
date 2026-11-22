## The skewness correction of the two-condition CR2 t (excess_skewness, skew_corrected_t in the core): a plain-R
## reference of kappa from the units' weights and residuals, Hall's transformation, and when it is (not) applied.

sk_rows <- function(m0, m1, seed, skew = TRUE, images = 2L) {
  set.seed(seed)
  m <- m0 + m1; unit <- rep(seq_len(m) - 1L, each = images); N <- length(unit)
  n <- round(runif(N, 20, 200)); v <- n * runif(N, 0.5, 2)
  eff <- if (skew) rexp(m)[unit + 1L] else rnorm(m)[unit + 1L]
  O <- n * (0.1 * eff + rnorm(N, 0, 0.05)); E <- rep(0, N)
  data.frame(img = seq_len(N) - 1L, O = O, E = E, n = n, v = v, unit = unit, group = as.integer(unit >= m0))
}

sk_hall <- function(t, k) t + k * t^2 / 3 + k^2 * t^3 / 27 + k / 6

## kappa from what the test returns: unit summaries (ybar = s / J), information J and tau2 give W = J / (1 + tau2 J)
sk_kappa_ref <- function(rows, u) {
  g <- tapply(rows$group, rows$unit, `[`, 1); id <- as.integer(names(g)) + 1L
  J <- u$unit_info[id]; ybar <- u$unit_summary[id]; W <- J / (1 + u$tau2 * J)
  S <- tapply(W, g, sum); beta <- tapply(W * ybar, g, sum) / S
  e <- ybar - beta[as.character(g)]; h <- W / S[as.character(g)]; z <- sqrt(W) * e
  w15 <- tapply(W^1.5, g, sum)
  unname(ref_shrunk_skewness(z, h, g) * (w15[2] / S[2]^3 - w15[1] / S[1]^3) / (1 / S[1] + 1 / S[2])^1.5)
}

test_that("kappa and the corrected p match the plain-R reference", {
  for (s in 1:4) {
    rows <- sk_rows(6 + 3 * s, 30, s)
    u <- spicyR:::stats_excess_test(rows, rows$unit, rows$group, max(rows$unit) + 1L, TRUE, "cr2", FALSE)
    expect_true(u$ok)
    k <- sk_kappa_ref(rows, u)
    expect_equal(u$skewness, k, tolerance = 1e-10)
    expect_equal(u$t, sk_hall(u$difference / u$se, k), tolerance = 1e-10)
    expect_equal(u$p, 2 * pt(-abs(sk_hall(u$difference / u$se, k)), u$df), tolerance = 1e-10)
    expect_equal(sum(u$influence^2), u$se^2, tolerance = 1e-10)   # se and influences are unchanged
  }
})

test_that("right-skewed unit effects with a small second condition give kappa > 0; swapping the conditions flips it", {
  rows <- sk_rows(40, 8, 11)
  m <- max(rows$unit) + 1L
  u <- spicyR:::stats_excess_test(rows, rows$unit, rows$group, m, TRUE, "cr2", FALSE)
  w <- spicyR:::stats_excess_test(rows, rows$unit, 1L - rows$group, m, TRUE, "cr2", FALSE)
  expect_gt(u$skewness, 0)
  expect_equal(w$skewness, -u$skewness, tolerance = 1e-12)
  expect_equal(w$p, u$p, tolerance = 1e-12)
})

test_that("equal conditions with equal weights: kappa is 0 and p is the uncorrected p", {
  m <- 10L; unit <- seq_len(2 * m) - 1L
  set.seed(3)
  rows <- data.frame(img = unit, O = rexp(2 * m), E = 0, n = 1, v = 1, unit = unit, group = as.integer(unit >= m))
  u <- spicyR:::stats_excess_test(rows, rows$unit, rows$group, 2L * m, FALSE, "cr2", FALSE)
  expect_identical(u$skewness, 0)
  expect_equal(u$p, 2 * pt(-abs(u$difference / u$se), u$df), tolerance = 1e-14)
})

test_that("Hartung-Knapp and the pooled dispersion are not corrected", {
  rows <- sk_rows(30, 8, 5)
  m <- max(rows$unit) + 1L
  for (u in list(spicyR:::stats_excess_test(rows, rows$unit, rows$group, m, TRUE, "hartung_knapp", FALSE),
                 spicyR:::stats_excess_test(rows, rows$unit, rows$group, m, TRUE, "cr2", TRUE))) {
    expect_identical(u$skewness, 0)
    expect_equal(u$p, 2 * pt(-abs(u$difference / u$se), u$df), tolerance = 1e-12)
  }
})

test_that("two units per condition give kappa 0 (the third moment has no residual df)", {
  rows <- data.frame(img = 0:3, O = c(1, 3, 2, 9), E = 0, n = 1, v = 1, unit = 0:3, group = c(0L, 0L, 1L, 1L))
  u <- spicyR:::stats_excess_test(rows, rows$unit, rows$group, 4L, FALSE, "cr2", FALSE)
  expect_identical(u$skewness, 0)
})

test_that("two units with nearly equal weights in each condition give kappa 0, not an unbounded one", {
  for (eps in c(0.1, 0.01, 0.001)) {
    rows <- data.frame(img = 0:3, O = c(0.3, -0.1, 0.05, -0.02) * c(1, 1 + eps, 9, 9 + eps), E = 0,
                       n = c(1, 1 + eps, 9, 9 + eps), v = c(1, 1 + eps, 9, 9 + eps), unit = 0:3, group = c(0L, 0L, 1L, 1L))
    u <- spicyR:::stats_excess_test(rows, rows$unit, rows$group, 4L, FALSE, "cr2", FALSE)
    expect_identical(u$skewness, 0)
    expect_equal(u$p, 2 * pt(-abs(u$difference / u$se), u$df), tolerance = 1e-12)
  }
})

test_that("Hall's transformation is monotone and odd in (t, kappa)", {
  t <- seq(-40, 40, by = 0.01)
  for (k in c(-1, -0.2, 0, 0.05, 0.3, 2)) {
    expect_true(all(diff(sk_hall(t, k)) >= 0))
    expect_equal(sk_hall(-t, -k), -sk_hall(t, k), tolerance = 1e-12)
  }
})

test_that("the skewness is shrunk by its own uncertainty: the correction fades in small studies", {
  # (that the package uses this shrinkage is checked against sk_kappa_ref above)
  set.seed(12)
  # normal residuals: the delta-method variance is about 6 / n
  z <- rnorm(1e5); V <- local({ s2 <- mean(z^2); t <- z / sqrt(s2); g <- mean(t^3); mean((t^3 - 3 * t - g - 1.5 * g * (t^2 - 1))^2) })
  expect_equal(V, 6, tolerance = 0.15)
  # the same skewed sample, replicated: the raw skewness is the same, its variance falls, so less shrinkage
  x <- rexp(12); f <- vapply(c(1, 4, 16, 256), function(r) { z <- rep(x - mean(x), r); h <- rep(1 / length(z), length(z))
    ref_shrunk_skewness(z, h) / ((sum(z^3) / sum((1 - h) * (1 - 2 * h))) / (sum(z^2) / sum(1 - h))^1.5) }, 0)
  expect_true(all(diff(f) > 0)); expect_gt(f[4], 0.95)
})
