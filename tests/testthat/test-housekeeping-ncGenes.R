test_that("Anticodon corresponds to a stop codon", {
    expect_error(
        tRNAGene(anticodon = 'UCA'),
        regexp = "The inserted anticodon corresponds to a stop codon"
    ) })


test_that("Anticodon doesn't match aminoacid", {
    expect_error(
        tRNAGene(anticodon = 'AAA', amino_acid = 'S'),
        regexp = "Anticodon and aminoacid are not consistent"
    ) })


test_that("Anticodon is not contained in mature tRNA sequence", {
    expect_error(
        tRNAGene(anticodon = 'UCA', tRNA_mature_sequence = 'UUUGCUA'),
        regexp = "Anticodon and mature tRNA sequence are not consistent"
    ) })


test_that("rRNA type is invalid", {
    expect_error(
        rRNAGene(rRNA_type = '4S'),
        regexp = "Please provide a valid human rRNA type"
    ) })


test_that("rRNA type is invalid w.r.t. chromosome", {
    expect_error(
        rRNAGene(rRNA_type = '18S', chromosome = 'chr7',
                strand = '+', gene_start = 20, gene_end = 40),
        regexp = "This rRNA type is encoded only by genes on chromosomes 13,"
    ) })



