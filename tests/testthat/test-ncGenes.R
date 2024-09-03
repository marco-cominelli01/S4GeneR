test_that("NonCodingGene has coding transcript among its
            alternative transcripts", {
    expect_error(
        {
            t1 <- tRNAGene()
            alternative_transcripts(t1, action = 'add') <- list(
                transcript_id = 'ENST12312312312', protein_coding = TRUE)
            },
        regexp = "A non-coding gene can't have a protein coding transcript"
    ) })
