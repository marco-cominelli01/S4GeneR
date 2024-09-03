test_that("miRNA's seed sequence is contained in mature miRNA sequence", {
    expect_error(
        miRNAGene(seed_sequence = 'AUAG', miRNA_mature_sequence = 'AAAAUAUA'),
        regexp = "The provided seed sequence is not compatible with"
    ) })
