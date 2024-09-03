test_that("Creation of valid empty objects via constructor is possible",
            expect_no_error(
                {
                    CodingGene()
                    tRNAGene()
                    rRNAGene()
                    miRNAGene()
                    siRNAGene()
                    piRNAGene()
                    snRNAGene()
                    snoRNAGene()
                    LongNonCodingGene()
                } ))


test_that("Creation of valid empty objects via new() is possible",
            expect_equal(
                {
                    x1 <- new("CodingGene")
                    x2 <- new("tRNAGene")
                    x3 <- new("rRNAGene")
                    x4 <- new("miRNAGene")
                    x5 <- new("siRNAGene")
                    x6 <- new("piRNAGene")
                    x7 <- new("snRNAGene")
                    x8 <- new("snoRNAGene")
                    x9 <- new("LongNonCodingGene")
                    obj_list <- list(x1,x2,x3,x4,x5,x6,x7,x8,x9)
                    all(unlist(lapply(obj_list, validObject)))
                }, expected = TRUE ))
