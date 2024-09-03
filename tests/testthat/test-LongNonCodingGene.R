test_that("lncRNA location class is valid", {
    expect_error(
        LongNonCodingGene(lncRNA_location_class = 'exogenous'),
        regexp = "Invalid lncRNA genomic location class"
    ) })


test_that("lncRNA functional class is valid", {
    expect_error(
        LongNonCodingGene(lncRNA_functional_class = 'transcription'),
        regexp = "Invalid lncRNA functional class"
    ) })
