## The inhomogeneous cross-K design (`sigma`). The reference numbers come from
## spicyR's spicyGLM(sigma = ) (gee branch) on the same seeded cells, so these
## pin the C++ core to the R implementation it ports.

make_compartment_cells <- function(seed = 3) {
  set.seed(seed)
  do.call(rbind, lapply(1:10, function(p) do.call(rbind, lapply(1:2, function(j) {
    n <- rpois(3, c(120, 160, 90)); tt <- rep(c("A", "B", "C"), n)
    left <- runif(length(tt)) < 0.75
    data.frame(x = ifelse(left, runif(length(tt), 0, 500), runif(length(tt), 500, 1000)),
               y = runif(length(tt), 0, 1000), cellType = tt, imageID = sprintf("i%02d_%d", p, j),
               subject = sprintf("s%02d", p), response = if (p <= 5) "alpha" else "beta",
               stringsAsFactors = FALSE)
  }))))
}

gee_reference <- list(
  convex = data.frame(
    from = c("A", "A", "A", "B", "B", "C"), to = c("A", "B", "C", "B", "C", "C"),
    coef_ref = c(-0.587922891303, -0.226469488559, -0.133166687398, -0.323012916693,
                 0.051282550928, -0.790426370318),
    coef_comp = c(-0.395933839417, -0.067043320317, 0.257108114667, -0.192219474991,
                  0.132259913738, -0.616378710901),
    p_value = c(0.460803743822, 0.027059375556, 0.237572746926, 0.290726701274,
                0.733921715634, 0.509681843156)),
  rectangle = data.frame(
    from = c("A", "A", "A", "B", "B", "C"), to = c("A", "B", "C", "B", "C", "C"),
    coef_ref = c(-0.538262029981, -0.200226575518, -0.064686045094, -0.280396931971,
                 0.109393579167, -0.775066374899),
    coef_comp = c(-0.356001742055, -0.023599133386, 0.290879931054, -0.179630687465,
                  0.182218261506, -0.597314466659),
    p_value = c(0.512930226224, 0.019193081727, 0.312366630399, 0.352409879406,
                0.761304909489, 0.514211731759)))

test_that("inhomogeneous fits reproduce spicyR's spicyGLM(sigma =)", {
  cells <- make_compartment_cells()
  for (window in names(gee_reference)) {
    ref <- gee_reference[[window]]
    out <- spicy_glm(cells, condition = "response", subject = "subject", r = 40, sigma = 150,
                     window = window)$results
    out <- out[match(paste(ref$from, ref$to), paste(out$from, out$to)), ]
    expect_equal(out$coef_ref, ref$coef_ref, tolerance = 1e-10)
    expect_equal(out$coef_comp, ref$coef_comp, tolerance = 1e-10)
    expect_equal(out$p_value, ref$p_value, tolerance = 1e-10)
  }
})

test_that("a disc covering every window, without edge correction, is the homogeneous fit", {
  cells <- make_compartment_cells()
  hom <- spicy_glm(cells, condition = "response", subject = "subject", r = 40)$results
  flat <- spicy_glm(cells, condition = "response", subject = "subject", r = 40, sigma = 1e7,
                    edge_correct = FALSE)$results
  ## cross pairs only: self-pairs of the homogeneous design count each cell as its own neighbour
  cross <- hom$from != hom$to
  flat <- flat[match(paste(hom$from, hom$to), paste(flat$from, flat$to)), ]
  expect_equal(flat$log_rate_ratio[cross], hom$log_rate_ratio[cross], tolerance = 1e-10)
  expect_equal(flat$p_value[cross], hom$p_value[cross], tolerance = 1e-8)
})

test_that("the inhomogeneous point estimate is direction-invariant", {
  cells <- make_compartment_cells()
  ab <- spicy_glm(cells, condition = "response", r = 40, sigma = 150, from = "A", to = "B")$results
  ba <- spicy_glm(cells, condition = "response", r = 40, sigma = 150, from = "B", to = "A")$results
  expect_equal(ab$log_rate_ratio, ba$log_rate_ratio, tolerance = 1e-12)
})

test_that("inhomogeneous diagnostics run, with leverage summing to one per group", {
  cells <- make_compartment_cells()
  out <- spicy_glm(cells, condition = "response", subject = "subject", r = 40, sigma = 150,
                   compute_diagnostics = TRUE)
  lev <- stats::aggregate(l_i ~ from + to + group, out$diagnostics$patient, sum)
  expect_equal(lev$l_i, rep(1, nrow(lev)), tolerance = 1e-12)
})

test_that("sigma is validated, and ignored for the binomial family", {
  cells <- make_compartment_cells()
  expect_error(spicy_glm(cells, condition = "response", r = 40, sigma = -1), "sigma")
  expect_error(spicy_glm(cells, condition = "response", r = 40, sigma = 150, min_lambda = 0),
               "min_lambda")
  expect_message(spicy_glm(cells, condition = "response", r = 40, sigma = 20, from = "A", to = "B"),
                 "no larger than")
  expect_message(b <- spicy_glm(cells, condition = "response", family = "binomial", k = 10,
                                sigma = 150, from = "A", to = "B"), "sigma")
  expect_equal(b$results, spicy_glm(cells, condition = "response", family = "binomial", k = 10,
                                    from = "A", to = "B")$results)
})
