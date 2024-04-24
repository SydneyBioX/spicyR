test_that(
    "the output is equal using segmented cells and SingleCellExperiment",
    {
        load(system.file("testdata/original_result.rda", package = "spicyR"))
        expect_equal(
            suppressWarnings(spicy(diabetesData_SCE,
                condition = "stage", subject = "case",
                from = "Tc", to = "Th"
            )),
            original_result
        )
    }
)
