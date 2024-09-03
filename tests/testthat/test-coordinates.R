test_that("Overlapping of exons coordinates", {
    expect_error(
        CodingGene(exons_starts = c(10,25), exons_ends = c(30,50),
                    chromosome = 'chr1', strand = '+'),
        regexp = "Exons can't overlap!"
    ) })


test_that("Valid chromosome", {
    expect_error(
        CodingGene(exons_starts = c(10,50), exons_ends = c(30,80),
                    chromosome = 'chrS', strand = '+'),
        regexp = "provide a valid human, or mitochondrial chromosome"
    ) })


test_that("Out of boundary mRNA coordinates w.r.t. gene", {
    expect_error(
        CodingGene(gene_start = 10, gene_end = 100,
                   chromosome = 'chrM', strand = '+',
                   mRNA_start = 30, mRNA_end = 120),
        regexp = "Out of boundaries coordinates"
    ) })


test_that("Out of boundary CDS coordinates w.r.t. gene", {
    expect_error(
        CodingGene(gene_start = 10, gene_end = 100,
                    chromosome = 'chrM', strand = '+',
                    cds_start = 30, cds_end = 120),
        regexp = "Out of boundaries coordinates"
    ) })


test_that("Out of boundary CDS coordinates w.r.t. mRNA", {
    expect_error(
        CodingGene(gene_start = 10, gene_end = 100,
                    chromosome = 'chrM', strand = '+',
                    mRNA_start = 20, mRNA_end = 80,
                    cds_start = 30, cds_end = 90),
        regexp = "Out of boundaries coordinates"
    ) })


test_that("CDS outside exons", {
    expect_error(
        CodingGene(gene_start = 10, gene_end = 100,
                    exons_starts = c(10,30,50), exons_ends = c(20, 40, 100),
                    chromosome = 'chrM', strand = '+',
                    mRNA_start = 20, mRNA_end = 80,
                    cds_start = 25, cds_end = 70),
        regexp = "must be inside the gene exons!"
    ) })


test_that("Gene starts and ends with exons", {
    expect_error(
        CodingGene(gene_start = 10, gene_end = 100,
                    exons_starts = c(11,30,50), exons_ends = c(20, 40, 100),
                    chromosome = 'chrM', strand = '+',
                    mRNA_start = 20, mRNA_end = 80,
                    cds_start = 25, cds_end = 70),
        regexp = "LAST exon must correspond to START and END"
    ) })

