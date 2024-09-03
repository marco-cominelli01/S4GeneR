#' @title Get Gene Product Length
#' @description Compute length of the gene product.
#' The specific attribute depends on the method implementation.
#' @param object Object of a class derived from `Gene`.
#' The specific type of object depends on the method implementation.
#' @return Gene product length.
#' @seealso
#' \link{lengthProduct,CodingGene-method},
#' \link{lengthProduct,tRNAGene-method},
#' \link{lengthProduct,tRNAGene-method},
#' \link{lengthProduct,miRNAGene-method},
#' \link{lengthProduct,siRNAGene-method},
#' \link{lengthProduct,piRNAGene-method},
#' \link{lengthProduct,snRNAGene-method},
#' \link{lengthProduct,snoRNAGene-method},
#' \link{lengthProduct,LongNonCodingGene-method},
#' @examples
#' # This is only the generic function, look in 'See Also' section
#' # for the method implementation.
#' # Example on a CodingGene object
#' cg1 <- CodingGene(protein_sequence = 'KRTA')
#' lengthProduct(cg1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setGeneric("lengthProduct", function(object) standardGeneric("lengthProduct"))


#' @title Get Protein Product Length
#' @description Compute length of protein_sequence attribute of `CodingGene`
#' class.
#' @param object Object of class `CodingGene`.
#' @return Protein product length.
#' @examples
#' cg1 <- CodingGene(protein_sequence = 'KRTA')
#' lengthProduct(cg1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("lengthProduct", "CodingGene", function(object) {
    if (length(object@protein_sequence) > 0) {
        output <- length(object@protein_sequence) }
    else {
        output <- NA }

    # It's possible to compute the protein length knowing exons and CDS
    # coordinates, by computing the intersection between them, summing the
    # widths of the resulting ranges and then dividing by 3. But, this
    # would make sense only if we know that the transcript includes entirely
    # (without any alternative splicing site event) the included exons.
    # Since alternative splicing site is not implemented, this constraint
    # is not considered.

    names(output) <- "encoded_protein_length"
    return(output) })


#' @title Get Mature tRNA Length
#' @description Compute length of tRNA_mature_sequence attribute of `tRNAGene`
#' class.
#' @param object Object of class `tRNAGene`.
#' @return Mature tRNA length.
#' @examples
#' t1 <- tRNAGene(tRNA_mature_sequence = 'UAUUAUAUAA')
#' lengthProduct(t1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("lengthProduct", "tRNAGene", function(object) {
    nonCodingLength(object, 'tRNA_mature_sequence',
                    'mature_tRNA_length') })


#' @title Get Mature rRNA Length
#' @description Compute length of rRNA_mature_sequence attribute of `rRNAGene`
#' class.
#' @param object Object of class `rRNAGene`.
#' @return Mature rRNA length.
#' @examples
#' r1 <- rRNAGene(rRNA_mature_sequence = 'UAUUAUAUAA')
#' lengthProduct(r1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("lengthProduct", "rRNAGene", function(object) {
    nonCodingLength(object, 'rRNA_mature_sequence',
                    'mature_rRNA_length') })


#' @title Get Mature miRNA Length
#' @description Compute length of miRNA_mature_sequence attribute of
#' `miRNAGene` class.
#' @param object Object of class `miRNAGene`.
#' @return Mature miRNA length.
#' @examples
#' mir1 <- miRNAGene(miRNA_mature_sequence = 'UAUUAUAUAA')
#' lengthProduct(mir1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("lengthProduct", "miRNAGene", function(object) {
    nonCodingLength(object, 'miRNA_mature_sequence',
                    'mature_miRNA_length') })


#' @title Get Mature siRNA Length
#' @description Compute length of siRNA_mature_sequence attribute of
#' `siRNAGene` class.
#' @param object Object of class `siRNAGene`.
#' @return Mature siRNA length.
#' @examples
#' sir1 <- siRNAGene(siRNA_mature_sequence = 'UAUUAUAUAA')
#' lengthProduct(sir1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("lengthProduct", "siRNAGene", function(object) {
    nonCodingLength(object, 'siRNA_mature_sequence',
                    'mature_siRNA_length') })


#' @title Get Mature piRNA Length
#' @description Compute length of piRNA_mature_sequence attribute of
#' `piRNAGene` class.
#' @param object Object of class `piRNAGene`.
#' @return Mature piRNA length.
#' @examples
#' pir1 <- piRNAGene(piRNA_mature_sequence = 'UAUUAUAUAA')
#' lengthProduct(pir1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("lengthProduct", "piRNAGene", function(object) {
    nonCodingLength(object, 'piRNA_mature_sequence',
                    'mature_piRNA_length') })


#' @title Get Mature snRNA Length
#' @description Compute length of snRNA_mature_sequence attribute of
#' `snRNAGene` class.
#' @param object Object of class `snRNAGene`.
#' @return Mature snRNA length.
#' @examples
#' sn1 <- snRNAGene(snRNA_mature_sequence = 'UAUUAUAUAA')
#' lengthProduct(sn1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("lengthProduct", "snRNAGene", function(object) {
    nonCodingLength(object, 'snRNA_mature_sequence',
                    'mature_snRNA_length') })


#' @title Get Mature snoRNA Length
#' @description Compute length of snoRNA_mature_sequence attribute of
#' `snoRNAGene` class.
#' @param object Object of class `snoRNAGene`.
#' @return Mature snoRNA length.
#' @examples
#' sno1 <- snoRNAGene(snoRNA_mature_sequence = 'UAUUAUAUAA')
#' lengthProduct(sno1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("lengthProduct", "snoRNAGene", function(object) {
    nonCodingLength(object, 'snoRNA_mature_sequence',
                    'mature_snoRNA_length') })


#' @title Get Mature lncRNA Length
#' @description Compute length of lncRNA_mature_sequence attribute of
#' `LongNonCodingGene` class.
#' @param object Object of class `LongNonCodingGene`.
#' @return Mature lncRNA length.
#' @examples
#' lnc1 <- LongNonCodingGene(lncRNA_mature_sequence = paste(rep('AUUUA', 50),
#' collapse = ''))
#' lengthProduct(lnc1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("lengthProduct", "LongNonCodingGene", function(object) {
    nonCodingLength(object, 'lncRNA_mature_sequence',
                    'mature_lncRNA_length') })
