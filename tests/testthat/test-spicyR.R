test_that(
    "the image method is unchanged from spicyR 1.x",
    {
        original_result <- readRDS(
            system.file("testdata/original_result.rds", package = "spicyR")
        )
        expect_equal(
            suppressWarnings(spicy(diabetesData,
                condition = "stage", subject = "case",
                from = "Tc", to = "Th", method = "image"
            )),
            original_result,
            tolerance = 0.01
        )
    }
)

test_that("the image method recycles a single from over several to", {
  res <- suppressWarnings(suppressMessages(spicy(diabetesData, condition = "stage", subject = "case",
                                                 from = "Tc", to = c("Th", "beta"), method = "image")))
  expect_equal(sort(rownames(res$p.value)), c("Tc__Th", "Tc__beta"))
})

test_that("the weight model fits a cell type paired with itself, as scam() does", {
    set.seed(3)
    x <- log10(rpois(60, 200) + 1)
    y <- 2 - 0.8 * x + rnorm(60, sd = 0.1)
    xnew <- seq(min(x), max(x), length.out = 25)
    z <- mpdWeightFit(y, x, x, xnew, xnew)
    expect_false(is.null(z))
    expect_true(all(diff(z) <= 1e-8))
    skip_if_not_installed("scam")
    fit <- scam::scam(y ~ s(x1, bs = "mpd") + s(x2, bs = "mpd"), data = data.frame(y = y, x1 = x, x2 = x))
    expected <- stats::predict(fit, data.frame(x1 = xnew, x2 = xnew))
    expect_equal(unname(z), as.vector(expected), tolerance = 1e-3)
})
