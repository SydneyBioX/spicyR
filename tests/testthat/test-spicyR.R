test_that(
    "the output is equal using segmented cells and SingleCellExperiment",
    {
        expect_equal(
            suppressWarnings(spicy(diabetesData_SCE,
                imageID = "imageID",
                condition = "stage", subject = "case",
                from = "Tc", to = "Th"
            )),
            suppressWarnings(spicy(diabetesData,
                condition = "stage", subject = "case",
                from = "Tc", to = "Th"
            ))
        )
    }
)
