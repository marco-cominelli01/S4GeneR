test_that("snRNA spliceosome complex validity", {
    expect_error(
        {
            snRNAGene(spliceosome_complex = 'U12')
        },
        regexp = "Please insert a valid type for the spliceosome complex"
    ) })
