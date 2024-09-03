test_that("RegulatoryNonCodingGene's target isn't a valid ENSEMBL ID", {
    expect_error(
        {
            mir1 <- miRNAGene(targets_ensembl_id = list('ENSG00012312312',
                                                        'ENST00012312313',
                                                        'ENSP00012312366'))
            targets_id(mir1, action = 'add') <- 'ENST0909090909A'
        },
        regexp = "The ID must start with ENSG/ENST/ENSP followed by exactly"
    ) })
