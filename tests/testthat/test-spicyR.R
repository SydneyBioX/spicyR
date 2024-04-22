# test_that(
#     "the output is equal using segmented cells and SingleCellExperiment",
#     {
#         expect_equal(
#             suppressWarnings(spicy(diabetesData_SCE,
#                 condition = "stage", subject = "case",
#                 from = "Tc", to = "Th"
#             )),
#             suppressWarnings(spicy(diabetesData,
#                 condition = "stage", subject = "case",
#                 from = "Tc", to = "Th"
#             ))
#         )
#     }
# )

test_that(
    "extractSpicyInfo makes an object identical to the original SegementedCells",
    {
        data(diabetesData_SCE)
        data(diabetesData)
        expect_equal(
            extractSpicyInfo(
                diabetesData_SCE,
                condition = "stage", subject = "case"
            ),
            diabetesData
        )
    }
)
