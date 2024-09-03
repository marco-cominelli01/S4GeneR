test_that("non-coding RNA mature length w.r.t. gene", {
    expect_error(
        tRNAGene(chromosome = 'chr1', strand = '+', gene_start = 20,
                gene_end = 25, tRNA_mature_sequence = 'AUACGGG'),
        regexp = "is not compatible with gene/exons/pre-ncRNA coordinates."
    ) })


test_that("non-coding RNA mature length w.r.t. exons", {
    expect_error(
        tRNAGene(chromosome = 'chr1', strand = '+', gene_start = 20,
                gene_end = 40, exons_starts = c(20,35), exons_ends = c(25,40),
                tRNA_mature_sequence = 'AUACGGGAUACGGG'),
        regexp = "is not compatible with gene/exons/pre-ncRNA coordinates."
    ) })


test_that("non-coding RNA mature length w.r.t. exons and pre-nc RNA", {
    expect_error(
        tRNAGene(chromosome = 'chr1', strand = '+', gene_start = 20,
                gene_end = 40, exons_starts = c(20,35), exons_ends = c(25,40),
                pre_ncRNA_start = 24, pre_ncRNA_end = 38,
                tRNA_mature_sequence = 'UCACCAC'),
        regexp = "is not compatible with gene/exons/pre-ncRNA coordinates."
    ) })


test_that("non-coding RNA mature length w.r.t. gene sequence", {
    expect_error(
        tRNAGene(gene_sequence = 'ATCGAA',
                tRNA_mature_sequence = 'UCACCAC'),
        regexp = "The latter can't be longer than the former."
    ) })


test_that("Small non-coding mature RNA longer than 200 nt", {
    expect_error(
        tRNAGene(chromosome = 'chr1', strand = '+',
                gene_start = 20, gene_end = 400,
                tRNA_mature_sequence = paste(rep('UCACCAC',40), collapse="")),
        regexp = "It should be shorter than 200 nucleotides."
    ) })


test_that("Long non-coding mature RNA shorter than 200 nt", {
    expect_error(
        LongNonCodingGene(chromosome = 'chr1', strand = '+',
                        gene_start = 20, gene_end = 400,
                        lncRNA_mature_sequence = paste(rep('UCACCAC',20),
                                                collapse="")),
        regexp = "It should be longer than 200 nucleotides."
    ) })


test_that("LongNonCodingGene's exons sum < 200 nt", {
    expect_error(
        LongNonCodingGene(chromosome = 'chr1', strand = '+',
                        gene_start = 20, gene_end = 400,
                        exons_starts = c(20, 350), exons_ends = c(90, 400)),
        regexp = "To encode a lncRNA, at least 200 nt are required."
    ) })


test_that("LongNonCodingGene's gene coordinates range < 200 nt", {
    expect_error(
        LongNonCodingGene(chromosome = 'chr1', strand = '+',
                        gene_start = 20, gene_end = 200),
        regexp = "To encode a lncRNA, at least 200 nt are required."
    ) })


test_that("LongNonCodingGene's gene sequence < 200 nt", {
    expect_error(
        LongNonCodingGene(gene_sequence = 'ATCATCA'),
        regexp = "A gene encoding for a lncRNA can't be shorter than 200 nt."
    ) })




