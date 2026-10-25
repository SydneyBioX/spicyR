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

test_that("the image method does not test a pair with the same value in every image", {
    # A fit to a constant response gives a p-value from rounding error alone, so the pair's p-value is NA (and it
    # stays out of the multiple-testing correction).
    set.seed(1)
    pheno <- data.frame(stage = rep(c("a", "b"), each = 12), case = rep(seq_len(8), each = 3))
    w <- runif(24, 0.5, 2)
    mem <- function(y) spatialMEM(y, cells = NULL, subject = "case", condition = "stage", covariates = NULL,
                                  weightFunction = w, cellCounts = NULL, pheno = pheno)
    expect_identical(mem(rep(-25, 24)), NA)
    expect_identical(mem(c(NA, rep(0, 23))), NA)
    expect_true(is.matrix(mem(rnorm(24))))
    lmFit <- function(y) spatialLM(y, from = 1, to = 1, cells = NULL, condition = "stage", covariates = NULL,
                                   weightFunction = w, cellCounts = matrix(1, 24, 1), pheno = pheno)
    expect_identical(lmFit(rep(-25, 24)), NA)
    expect_s3_class(lmFit(rnorm(24)), "lm")
    p <- cleanMEM(list(a__b = mem(rep(-25, 24)), b__a = mem(rnorm(24))))$p.value
    expect_true(is.na(p["a__b", "conditionb"]) && !is.na(p["b__a", "conditionb"]))
})

test_that("the weight model follows scam()'s smoothing-parameter search", {
    # The GCV surface is flat near its minimum, so a search that strays from
    # scam()'s path (e.g. with a finite-difference gradient) stops elsewhere
    # within tolerance and moves the weights; the fit must match scam() closely.
    skip_if_not_installed("scam")
    res <- suppressWarnings(suppressMessages(spicy(diabetesData,
        condition = "stage", subject = "case", method = "image"
    )))
    pa <- do.call(cbind, res$pairwiseAssoc)
    rS <- as.vector(apply(pa, 2, function(x) {
        if (sum(!is.na(x)) > 1 && stats::sd(x, na.rm = TRUE) > 0) (x - mean(x, na.rm = TRUE))^2 else rep(NA, length(x))
    }))
    ok <- !is.na(rS)
    y <- log10(rS[ok] + 1)
    x1 <- log10(as.vector(res$nCells[, res$comparisons$from])[ok] + 1)
    x2 <- log10(as.vector(res$nCells[, res$comparisons$to])[ok] + 1)
    z <- mpdWeightFit(y, x1, x2, x1, x2)
    fit <- scam::scam(y ~ s(x1, bs = "mpd") + s(x2, bs = "mpd"))
    expect_equal(unname(z), unname(stats::fitted(fit)), tolerance = 1e-10)
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
