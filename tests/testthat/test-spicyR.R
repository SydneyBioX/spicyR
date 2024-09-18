test_that(
    "the output is equal to previous version",
    {
        original_result <- readRDS(
            system.file("testdata/original_result.rds", package = "spicyR")
        )
        expect_equal(
            suppressWarnings(spicy(diabetesData,
                condition = "stage", subject = "case",
                from = "Tc", to = "Th"
            )),
            original_result
        )
    }
)
