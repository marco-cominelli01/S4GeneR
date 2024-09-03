test_that("snoRNA box type is valid", {
    expect_error(
        snoRNAGene(snoRNA_box_type = 'ABCD'),
        regexp = "Invalid snoRNA box type"
    ) })
