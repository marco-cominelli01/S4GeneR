test_that("Gene sequence length compatible with gene and exons coordinates", {
    expect_error(
        CodingGene(chromosome = 'chr1', strand = '+',
                    gene_start = 20, gene_end = 30,
                    gene_sequence = 'ATCGATTAAGGG'),
        regexp = "Length of gene sequence and gene/exons coordinates are not"
    )
    expect_no_error(
        CodingGene(chromosome = 'chr1', strand = '+',
                    gene_start = 20, gene_end = 30,
                    gene_sequence = 'ATCGATTAAGG')
    )
    expect_error(
        CodingGene(chromosome = 'chr1', strand = '+',
                    exons_starts = 20, exons_ends = 30,
                    gene_sequence = 'ATCGATTAAGGG'),
        regexp = "Length of gene sequence and gene/exons coordinates are not"
    )
    expect_no_error(
        CodingGene(chromosome = 'chr1', strand = '+',
                    exons_starts = 20, exons_ends = 30,
                    gene_sequence = 'ATCGATTAAGG'),
    )})


test_that("Adding already present alternative transcript",
    expect_error(
        {
            cg1 <- CodingGene()
            alternative_transcripts(cg1, action = 'add') <- list(
                transcript_id = 'ENST12312312312', protein_coding = TRUE)
            alternative_transcripts(cg1, action = 'add') <- list(
                transcript_id = 'ENST12312312312', protein_coding = FALSE)
        },
        regexp = "Transcript already present among the alternative"
    ))
