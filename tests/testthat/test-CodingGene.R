test_that("Protein length validity w.r.t. gene sequence", {
    expect_error(
        CodingGene(protein_sequence = 'KRTK',
                    gene_sequence = 'AGGCCCGGGCCG'),
        regexp = "Any protein coding gene"
    ) })


test_that("Protein length validity w.r.t. CDS and exons coordinates", {
    expect_error(
        CodingGene(chromosome = 'chr1', strand = '+',
                    exons_starts = c(10,40,50), exons_ends = c(20,45,55),
                    cds_start = 15, cds_end = 43,
                    protein_sequence = 'KRTK'),
        regexp = "intersection between CDS and exons, divided by 3"
    )
    expect_no_error(
        CodingGene(chromosome = 'chr1', strand = '+',
                    exons_starts = c(10,40,50), exons_ends = c(20,45,55),
                    cds_start = 15, cds_end = 43,
                    protein_sequence = 'KRK')
    )})
