test_that("Product length of CodingGene object", {
    expected <- 4
    names(expected) <- 'encoded_protein_length'
    expect_equal(
        lengthProduct(CodingGene(protein_sequence = 'KRTK')),
        expected
    ) })


test_that("Product length of tRNAGene object", {
    expected <- 5
    names(expected) <- 'mature_tRNA_length'
    expect_equal(
        lengthProduct(tRNAGene(tRNA_mature_sequence = 'AUCAG')),
        expected
    ) })


test_that("Product length of rRNAGene object", {
    expected <- 7
    names(expected) <- 'mature_rRNA_length'
    expect_equal(
        lengthProduct(rRNAGene(rRNA_mature_sequence = 'AGGAGAG')),
        expected
    ) })


test_that("Product length of miRNAGene object", {
    expected <- 4
    names(expected) <- 'mature_miRNA_length'
    expect_equal(
        lengthProduct(miRNAGene(miRNA_mature_sequence = 'AUUA')),
        expected
    ) })


test_that("Product length of siRNAGene object", {
    expected <- 5
    names(expected) <- 'mature_siRNA_length'
    expect_equal(
        lengthProduct(siRNAGene(siRNA_mature_sequence = 'UUAUA')),
        expected
    ) })


test_that("Product length of piRNAGene object", {
    expected <- 10
    names(expected) <- 'mature_piRNA_length'
    expect_equal(
        lengthProduct(piRNAGene(piRNA_mature_sequence = 'AAAAAAAAAA')),
        expected
    ) })


test_that("Product length of snRNAGene object", {
    expected <- 1
    names(expected) <- 'mature_snRNA_length'
    expect_equal(
        lengthProduct(snRNAGene(snRNA_mature_sequence = 'G')),
        expected
    ) })


test_that("Product length of snoRNAGene object", {
    expected <- NA
    names(expected) <- 'mature_snoRNA_length'
    expect_equal(
        lengthProduct(snoRNAGene()),
        expected
    ) })


test_that("Product length of LongNonCodingGene object", {
    expected <- 300
    names(expected) <- 'mature_lncRNA_length'
    expect_equal(
        lengthProduct(LongNonCodingGene(lncRNA_mature_sequence =
                                paste(rep("AGAGA",60), collapse = ""))),
        expected
    ) })
