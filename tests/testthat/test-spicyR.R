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
