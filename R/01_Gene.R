#' @title Gene Virtual Class
#' @description This virtual class is the one from which all the others
#' derives.
#' @slot gene_ensembl_id `character`. Gene ENSEMBL ID.
#' @slot hugo_symbol `character`. Gene Hugo symbol.
#' @slot gene_complete_name `character`. Gene complete name.
#' @slot gene_description `character`. Gene description.
#' @slot gene_coordinates `GRanges`. Gene genomic coordinates.
#' @slot gene_sequence `DNAString`. Gene sequence.
#' @slot exons_coordinates `GRanges`. Exons genomic coordinates.
#' @slot alternative_transcripts `data.frame`.
#' Alternative transcripts.
#' @return No objects can be created directly from this class since
#' it's virtual.
#' @seealso
#' \linkS4class{CodingGene}, \linkS4class{NonCodingGene}
#' @examples
#' # The Gene class is virtual, so no objects can be created directly
#' # from this class. It needs to be extended by other classes.
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setClass("Gene",
        contains = "VIRTUAL",
        slots = c(
            gene_ensembl_id = "character",
            hugo_symbol = "character",
            gene_complete_name = "character",
            gene_description = "character",
            gene_coordinates = "GRanges",
            gene_sequence = "DNAString",
            exons_coordinates = "GRanges",
            alternative_transcripts = "data.frame"
        ),
        prototype = list(
            gene_ensembl_id = "ENSG00000000000",
            hugo_symbol = NA_character_,
            gene_complete_name = NA_character_,
            gene_description = NA_character_,
            gene_coordinates = GRanges(),
            gene_sequence = DNAString(),
            exons_coordinates = GRanges(),
            alternative_transcripts = data.frame(
                transcript_ensembl_id = 'ENST00000000000',
                protein_coding = FALSE)
        ))


#' @title Validity Check for Gene Virtual Class
#' @description The function checks that: each slot, except the special ones,
#' is of length 1 (or empty); each ID is a valid ENSEMBL ID; all coordinates
#' are consistent w.r.t. each other and that gene sequence is consistent
#' with the provided coordinates.
#' @param object Object of class `Gene` (or derived classes).
#' @return `TRUE` if the object is valid, otherwise the function
#' throws an error via `stop()`.
#' @seealso
#' \linkS4class{Gene}
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @name Gene-validity

setValidity("Gene", function(object){
    ### 1. Check on slots length (they can contain at most 1 element, except
    ###    the special ones).
    special_slots <- c( "exons_coordinates","alternative_transcripts",
                        "gene_sequence" )
    check_slots_length( obj = object, class_obj = "Gene",
                        mother_class_obj = NA, special_slots = special_slots )
    ### 2. Check the validity of ENSEMBL IDs.
    ### a) Gene ID
    gene_id <- object@gene_ensembl_id
    check_ensembl_id(id = gene_id, prefix = 'ENSG', type = 'gene')
    ### b) Transcript ID inside the slot 'additional_alternative_transcript'.
    msg <- " in 'alternative_transcripts' slot"
    alt_transcript_ids <- object@alternative_transcripts$transcript_ensembl_id
    alt_transcript_ids <- as.list(alt_transcript_ids)
    lapply(alt_transcript_ids, check_ensembl_id, prefix = 'ENST',
            type = 'transcript', add_msg = msg)
    ### 3. Check that gene and exons coordinates (if defined) share both
    ###    chromosome and strand.
    check_coordinates_consistency(granges_list = list(
        gene_coordinates = object@gene_coordinates,
        exons_coordinates = object@exons_coordinates ))
    ### 4. Check that exons coordinates are valid (distinct exons should not
    ### overlap) and that exons make sense w.r.t gene coordinates.
    check_exons_coordinates(exons = object@exons_coordinates,
                            gene = object@gene_coordinates)
    ### 5. Check that length of gene sequence is compatible with coordinates.
    ###    If gene coordinates are defined, gene sequence must be of the same
    ###    length w.r.t the width of the gene coordinates range. Same holds
    ###    if only exons coordinates are defined.
    coordinates <- c(object@gene_coordinates,
                    object@exons_coordinates)
    # If at least one of the two is not empty.
    if (!isEmpty(coordinates)) {
        exons_width <- width(range(object@exons_coordinates))
        gene_width <- width(object@gene_coordinates)
        gene_length <- max(exons_width, gene_width)
        gene_sequence_length <- length(object@gene_sequence)
        if (gene_sequence_length && gene_sequence_length != gene_length) {
            err_msg <- paste("Length of gene sequence and gene/exons",
                        "coordinates are not coherent: the length of the",
                        "sequence should be (gene_end - gene_start + 1) or,",
                        "if just exons coordinates are defined,",
                        "(last_exon_end - first_exon_start + 1).")
            stop(err_msg)  } }
    return(TRUE) })


# GETTERs for Gene

#' @title Get Gene ENSEMBL ID
#' @description Getter for gene_ensembl_id attribute of `Gene` class
#' (and derived classes).
#' @param object Object of class `Gene` (and derived classes).
#' @return Value of gene_ensembl_id.
#' @seealso
#' \link{gene_id,Gene-method}
#' @examples
#' # This is only the generic function, look in 'See Also' section
#' # for the method implementation.
#' cg1 <- CodingGene(gene_ensembl_id = 'ENSG12312312300')
#' gene_id(cg1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setGeneric("gene_id", function(object) standardGeneric("gene_id") )


#' @title Get Gene ENSEMBL ID
#' @description Getter for gene_ensembl_id attribute of `Gene` class
#' (and derived classes).
#' @param object Object of class `Gene` (and derived classes).
#' @return Value of gene_ensembl_id.
#' @examples
#' # Since Gene is a virtual class, the example is made using a derived
#' # class.
#' cg1 <- CodingGene(gene_ensembl_id = 'ENSG12312312300')
#' gene_id(cg1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("gene_id", "Gene", function(object) { object@gene_ensembl_id })


#' @title Get Gene Hugo symbol
#' @description Getter for hugo_symbol attribute of `Gene` class
#' (and derived classes).
#' @param object Object of class `Gene` (and derived classes).
#' @return Value of hugo_symbol.
#' @seealso
#' \link{hugo_symbol,Gene-method}
#' @examples
#' # This is only the generic function, look in 'See Also' section
#' # for the method implementation.
#' cg1 <- CodingGene(hugo_symbol = 'TP53')
#' hugo_symbol(cg1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setGeneric("hugo_symbol", function(object) standardGeneric("hugo_symbol") )


#' @title Get Gene Hugo symbol
#' @description Getter for hugo_symbol attribute of `Gene` class
#' (and derived classes).
#' @param object Object of class `Gene` (and derived classes).
#' @return Value of hugo_symbol.
#' @examples
#' # Since Gene is a virtual class, the example is made using a derived
#' # class.
#' cg1 <- CodingGene(hugo_symbol = 'TP53')
#' hugo_symbol(cg1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("hugo_symbol", "Gene", function(object) { object@hugo_symbol })


#' @title Get Gene Complete Name
#' @description Getter for gene_complete_name attribute of `Gene` class
#' (and derived classes).
#' @param object Object of class `Gene` (and derived classes).
#' @return Value of gene_complete_name.
#' @seealso
#' \link{gene_complete_name,Gene-method}
#' @examples
#' # This is only the generic function, look in 'See Also' section
#' # for the method implementation.
#' cg1 <- CodingGene(gene_complete_name = 'Tumor Protein 53')
#' gene_complete_name(cg1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setGeneric("gene_complete_name",
            function(object) standardGeneric("gene_complete_name") )


#' @title Get Gene Complete Name
#' @description Getter for gene_complete_name attribute of `Gene` class
#' (and derived classes).
#' @param object Object of class `Gene` (and derived classes).
#' @return Value of gene_complete_name.
#' @examples
#' # Since Gene is a virtual class, the example is made using a derived
#' # class.
#' cg1 <- CodingGene(gene_complete_name = 'Tumor Protein 53')
#' gene_complete_name(cg1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("gene_complete_name", "Gene", function(object) {
    object@gene_complete_name })


#' @title Get Gene Description
#' @description Getter for gene_description attribute of `Gene` class
#' (and derived classes).
#' @param object Object of class `Gene` (and derived classes).
#' @return Value of gene_description.
#' @seealso
#' \link{gene_description,Gene-method}
#' @examples
#' # This is only the generic function, look in 'See Also' section
#' # for the method implementation.
#' cg1 <- CodingGene(gene_description = 'TP53 is a tumor suppressor.')
#' gene_description(cg1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setGeneric("gene_description",
            function(object) standardGeneric("gene_description") )


#' @title Get Gene Description
#' @description Getter for gene_description attribute of `Gene` class
#' (and derived classes).
#' @param object Object of class `Gene` (and derived classes).
#' @return Value of gene_description.
#' @examples
#' # Since Gene is a virtual class, the example is made using a derived
#' # class.
#' cg1 <- CodingGene(gene_description = 'TP53 is a tumor suppressor.')
#' gene_description(cg1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("gene_description", "Gene", function(object) {
    object@gene_description })


#' @title Get Gene Genomic Coordinates
#' @description Getter for gene_coordinates attribute of `Gene` class
#' (and derived classes).
#' @param object Object of class `Gene` (and derived classes).
#' @return `GRanges` object representing gene coordinates.
#' @seealso
#' \link{gene_coordinates,Gene-method}
#' @examples
#' # This is only the generic function, look in 'See Also' section
#' # for the method implementation.
#' cg1 <- CodingGene(chromosome = 'chr1', strand = '+', gene_start = 210,
#' gene_end = 800)
#' gene_coordinates(cg1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setGeneric("gene_coordinates",
            function(object) standardGeneric("gene_coordinates") )


#' @title Get Gene Genomic Coordinates
#' @description Getter for gene_coordinates attribute of `Gene` class
#' (and derived classes).
#' @param object Object of class `Gene` (and derived classes).
#' @return `GRanges` object representing gene coordinates.
#' @examples
#' # Since Gene is a virtual class, the example is made using a derived
#' # class.
#' cg1 <- CodingGene(chromosome = 'chr1', strand = '+', gene_start = 210,
#' gene_end = 800)
#' gene_coordinates(cg1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("gene_coordinates", "Gene", function(object) {
    object@gene_coordinates })


#' @title Get Gene Sequence
#' @description Getter for gene_sequence attribute of `Gene` class
#' (and derived classes).
#' @param object Object of class `Gene` (and derived classes).
#' @return Gene sequence.
#' @seealso
#' \link{gene_sequence,Gene-method}
#' @examples
#' # This is only the generic function, look in 'See Also' section
#' # for the method implementation.
#' cg1 <- CodingGene(gene_sequence = 'ATCATCATCACCATC')
#' gene_sequence(cg1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setGeneric("gene_sequence",
            function(object) standardGeneric("gene_sequence") )


#' @title Get Gene Sequence
#' @description Getter for gene_sequence attribute of `Gene` class
#' (and derived classes).
#' @param object Object of class `Gene` (and derived classes).
#' @return Gene sequence.
#' @examples
#' # Since Gene is a virtual class, the example is made using a derived
#' # class.
#' cg1 <- CodingGene(gene_sequence = 'ATCATCATCACCATC')
#' gene_sequence(cg1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("gene_sequence", "Gene", function(object) {
    object@gene_sequence })


#' @title Get Exons Genomic Coordinates
#' @description Getter for exons_coordinates attribute of `Gene` class
#' (and derived classes).
#' @param object Object of class `Gene` (and derived classes).
#' @return `GRanges` object representing exons coordinates.
#' @seealso
#' \link{exons_coordinates,Gene-method}
#' @examples
#' # This is only the generic function, look in 'See Also' section
#' # for the method implementation.
#' cg1 <- CodingGene(chromosome = 'chr1', strand = '+',
#' exons_starts = c(10,50), exons_ends = c(30, 90))
#' exons_coordinates(cg1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setGeneric("exons_coordinates",
            function(object) standardGeneric("exons_coordinates") )


#' @title Get Exons Genomic Coordinates
#' @description Getter for exons_coordinates attribute of `Gene` class
#' (and derived classes).
#' @param object Object of class `Gene` (and derived classes).
#' @return `GRanges` object representing exons coordinates.
#' @examples
#' # Since Gene is a virtual class, the example is made using a derived
#' # class.
#' cg1 <- CodingGene(chromosome = 'chr1', strand = '+',
#' exons_starts = c(10,50), exons_ends = c(30, 90))
#' exons_coordinates(cg1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("exons_coordinates", "Gene", function(object) {
    object@exons_coordinates })


#' @title Get Alternative Transcripts
#' @description Getter for alternative_transcripts attribute of `Gene` class
#' (and derived classes). This attribute can't be modified at construction
#' time because transcripts can be added/removed only after the creation of
#' the object by using the dedicated setter.
#' @param object Object of class `Gene` (and derived classes).
#' @return Dataframe with transcript_id and protein_coding (`logical`)
#' columns.
#' @seealso
#' \link{alternative_transcripts,Gene-method}
#' @examples
#' # This is only the generic function, look in 'See Also' section
#' # for the method implementation.
#' cg1 <- CodingGene()
#' alternative_transcripts(cg1)
#' # The returned row is just a placeholder (it's removed after the first
#' # modification done by the user).
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setGeneric("alternative_transcripts",
            function(object) standardGeneric("alternative_transcripts") )


#' @title Get alternative transcripts dataframe
#' @description Getter for alternative_transcripts attribute of `Gene` class
#' (and derived classes). This attribute can't be modified at construction
#' time but transcripts can be added/removed only after the creation of the
#' object by using the dedicated setter.
#' @param object Object of class `Gene` (and derived classes).
#' @return Dataframe with transcript_id and protein_coding (`logical`)
#' columns.
#' @examples
#' # Since Gene is a virtual class, the example is made using a derived
#' # class.
#' cg1 <- CodingGene()
#' alternative_transcripts(cg1)
#' # The returned row is just a placeholder (it's removed after the first
#' # modification done by the user).
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("alternative_transcripts", "Gene", function(object) {
    object@alternative_transcripts })


# SETTERs for Gene


#' @title Set Gene ENSEMBL ID
#' @description Setter for gene_ensembl_id attribute of `Gene` class
#' (and derived classes).
#' @param object Object of class `Gene` (and derived classes).
#' @param value New value for gene_ensembl_id.
#' @return The modified object.
#' @seealso
#' \link{gene_id<-,Gene-method}
#' @examples
#' # This is only the generic function, look in 'See Also' section
#' # for the method implementation.
#' cg1 <- CodingGene(gene_ensembl_id = 'ENSG12312312300')
#' gene_id(cg1)
#' gene_id(cg1) <- 'ENSG00099988812'
#' gene_id(cg1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setGeneric("gene_id<-", function(object, value) standardGeneric("gene_id<-") )


#' @title Set Gene ENSEMBL ID
#' @description Setter for gene_ensembl_id attribute of `Gene` class
#'  (and derived classes).
#' @param object Object of class `Gene` (and derived classes).
#' @param value New value for gene_ensembl_id.
#' @return The modified object.
#' @examples
#' # Since Gene is a virtual class, the example is made using a derived
#' # class.
#' cg1 <- CodingGene(gene_ensembl_id = 'ENSG12312312300')
#' gene_id(cg1)
#' gene_id(cg1) <- 'ENSG00099988812'
#' gene_id(cg1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("gene_id<-", "Gene", function(object, value) {
    object@gene_ensembl_id <- value
    validObject(object)
    return(object) })


#' @title Set Gene Hugo symbol
#' @description Setter for hugo_symbol attribute of `Gene` class
#' (and derived classes).
#' @param object Object of class `Gene` (and derived classes).
#' @param value New value for hugo_symbol.
#' @return The modified object.
#' @seealso
#' \link{hugo_symbol<-,Gene-method}
#' @examples
#' # This is only the generic function, look in 'See Also' section
#' # for the method implementation.
#' cg1 <- CodingGene(hugo_symbol = 'TP53')
#' hugo_symbol(cg1)
#' hugo_symbol(cg1) <- 'BRCA1'
#' hugo_symbol(cg1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setGeneric("hugo_symbol<-",
            function(object, value) standardGeneric("hugo_symbol<-") )


#' @title Set Gene Hugo symbol
#' @description Setter for hugo_symbol attribute of `Gene` class
#' (and derived classes).
#' @param object Object of class `Gene` (and derived classes).
#' @param value New value for hugo_symbol.
#' @return The modified object.
#' @examples
#' # Since Gene is a virtual class, the example is made using a derived
#' # class.
#' cg1 <- CodingGene(hugo_symbol = 'TP53')
#' hugo_symbol(cg1)
#' hugo_symbol(cg1) <- 'BRCA1'
#' hugo_symbol(cg1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("hugo_symbol<-", "Gene", function(object, value) {
    object@hugo_symbol <- value
    validObject(object)
    return(object) })


#' @title Set Gene Complete Name
#' @description Setter for gene_complete_name attribute of `Gene` class
#' (and derived classes).
#' @param object Object of class `Gene` (and derived classes).
#' @param value New value for gene_complete_name.
#' @return The modified object.
#' @seealso
#' \link{gene_complete_name<-,Gene-method}
#' @examples
#' # This is only the generic function, look in 'See Also' section
#' # for the method implementation.
#' cg1 <- CodingGene(gene_complete_name = 'Tumor Protein 53')
#' gene_complete_name(cg1)
#' gene_complete_name(cg1) <- 'Breast cancer type 1 susceptibility protein'
#' gene_complete_name(cg1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setGeneric("gene_complete_name<-",
            function(object, value) standardGeneric("gene_complete_name<-") )


#' @title Set Gene Complete Name
#' @description Setter for gene_complete_name attribute of `Gene` class
#' (and derived classes).
#' @param object Object of class `Gene` (and derived classes).
#' @param value New value for gene_complete_name.
#' @return The modified object.
#' @examples
#' # Since Gene is a virtual class, the example is made using a derived
#' # class.
#' cg1 <- CodingGene(gene_complete_name = 'Tumor Protein 53')
#' gene_complete_name(cg1)
#' gene_complete_name(cg1) <- 'Breast cancer type 1 susceptibility protein'
#' gene_complete_name(cg1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("gene_complete_name<-", "Gene", function(object, value) {
    object@gene_complete_name <- value
    validObject(object)
    return(object) })


#' @title Set Gene Description
#' @description Setter for gene_description attribute of `Gene` class
#' (and derived classes).
#' @param object Object of class `Gene` (and derived classes).
#' @param value New value for gene_description.
#' @return The modified object.
#' @seealso
#' \link{gene_description<-,Gene-method}
#' @examples
#' # This is only the generic function, look in 'See Also' section
#' # for the method implementation.
#' cg1 <- CodingGene(gene_description = 'TP53 is a tumor suppressor gene')
#' gene_description(cg1)
#' gene_description(cg1) <- 'BRCA1 is involved in breast cancer onset'
#' gene_description(cg1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setGeneric("gene_description<-",
            function(object, value) standardGeneric("gene_description<-") )

#' @title Set Gene Description
#' @description Setter for gene_description attribute of `Gene` class
#' (and derived classes).
#' @param object Object of class `Gene` (and derived classes).
#' @param value New value for gene_description.
#' @return The modified object.
#' @examples
#' # Since Gene is a virtual class, the example is made using a derived
#' # class.
#' cg1 <- CodingGene(gene_description = 'TP53 is a tumor suppressor gene')
#' gene_description(cg1)
#' gene_description(cg1) <- 'BRCA1 is involved in breast cancer onset'
#' gene_description(cg1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("gene_description<-", "Gene", function(object, value) {
    object@gene_description <- value
    validObject(object)
    return(object) })

#' @title Set Gene Genomic Coordinates
#' @description Setter for gene_coordinates attribute of `Gene` class
#' (and derived classes).
#' @param object Object of class `Gene` (and derived classes).
#' @param value Named list with chromosome, strand, start and end values
#' to build the new gene coordinates.
#' @return The modified object.
#' @seealso
#' \link{gene_coordinates<-,Gene-method}
#' @examples
#' # This is only the generic function, look in 'See Also' section
#' # for the method implementation.
#' cg1 <- CodingGene(chromosome = 'chr1', strand = '+', gene_start = 20,
#' gene_end = 300)
#' gene_coordinates(cg1)
#' gene_coordinates(cg1) <- list(chromosome = 'chr2', strand = '+',
#' start = 30, end = 50)
#' gene_coordinates(cg1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setGeneric("gene_coordinates<-",
            function(object, value) standardGeneric("gene_coordinates<-") )


#' @title Set Gene Genomic Coordinates
#' @description Setter for gene_coordinates attribute of `Gene` class
#' (and derived classes).
#' @param object Object of class `Gene` (and derived classes).
#' @param value Named list with chromosome, strand, start and end values
#' to build the new gene coordinates.
#' @return The modified object.
#' @examples
#' # Since Gene is a virtual class, the example is made using a derived
#' # class.
#' cg1 <- CodingGene(chromosome = 'chr1', strand = '+', gene_start = 20,
#' gene_end = 300)
#' gene_coordinates(cg1)
#' gene_coordinates(cg1) <- list(chromosome = 'chr2', strand = '+',
#' start = 30, end = 50)
#' gene_coordinates(cg1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("gene_coordinates<-", "Gene", function(object, value) {
    list_names <- c("chromosome", "strand", "start", "end")
    msg <- paste0("e.g. 'gene_coordinates(g1) <- list(chromosome = 'chr1',",
                    " strand = '+', start = 10, end = 260)'.")
    object@gene_coordinates <- create_new_coordinates(value = value,
                                                    list_names = list_names,
                                                    msg = msg)
    validObject(object)
    return(object) })


#' @title Set Gene Sequence
#' @description Setter for gene_sequence attribute of `Gene` class
#' (and derived classes).
#' @param object Object of class `Gene` (and derived classes).
#' @param value New value for gene_sequence.
#' @return The modified object.
#' @seealso
#' \link{gene_sequence<-,Gene-method}
#' @examples
#' # This is only the generic function, look in 'See Also' section
#' # for the method implementation.
#' cg1 <- CodingGene(gene_sequence = 'ATCCAT')
#' gene_sequence(cg1)
#' gene_sequence(cg1) <- 'AAAAAAAA'
#' gene_sequence(cg1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setGeneric("gene_sequence<-",
            function(object, value) standardGeneric("gene_sequence<-") )


#' @title Set Gene Sequence
#' @description Setter for gene_sequence attribute of `Gene` class
#' (and derived classes).
#' @param object Object of class `Gene` (and derived classes).
#' @param value New value for gene_sequence.
#' @return The modified object.
#' @examples
#' # Since Gene is a virtual class, the example is made using a derived
#' # class.
#' cg1 <- CodingGene(gene_sequence = 'ATCCAT')
#' gene_sequence(cg1)
#' gene_sequence(cg1) <- 'AAAAAAAA'
#' gene_sequence(cg1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("gene_sequence<-", "Gene", function(object, value) {
    object@gene_sequence <- DNAString(as.character(value))
    validObject(object)
    return(object) })


#' @title Set Exons Genomic Coordinates
#' @description Setter for exons_coordinates attribute of `Gene` class
#' (and derived classes).
#' @param object Object of class `Gene` (and derived classes).
#' @param value Named list with chromosome, strand, starts and ends values
#' to build the new exons coordinates.
#' @return The modified object.
#' @seealso
#' \link{exons_coordinates<-,Gene-method}
#' @examples
#' # This is only the generic function, look in 'See Also' section
#' # for the method implementation.
#' cg1 <- CodingGene(chromosome = 'chr1', strand = '+',
#' exons_starts = c(20,60), exons_ends = c(30,70))
#' exons_coordinates(cg1)
#' exons_coordinates(cg1) <- list(chromosome = 'chr2', strand = '+',
#' starts = c(90,150), ends = c(100,400))
#' exons_coordinates(cg1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setGeneric("exons_coordinates<-",
            function(object, value) standardGeneric("exons_coordinates<-"))

#' @title Set Exons Genomic Coordinates
#' @description Setter for exons_coordinates attribute of `Gene` class
#' (and derived classes).
#' @param object Object of class `Gene` (and derived classes).
#' @param value Named list with chromosome, strand, starts and ends values
#' to build the new exons coordinates.
#' @return The modified object.
#' @examples
#' # Since Gene is a virtual class, the example is made using a derived
#' # class.
#' cg1 <- CodingGene(chromosome = 'chr1', strand = '+',
#' exons_starts = c(20,60), exons_ends = c(30,70))
#' exons_coordinates(cg1)
#' exons_coordinates(cg1) <- list(chromosome = 'chr2', strand = '+',
#' starts = c(90,150), ends = c(100,400))
#' exons_coordinates(cg1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("exons_coordinates<-", "Gene", function(object, value) {
    if (length(value$starts) != length(value$ends))
        stop("Each exon start must match an exon end.")

    list_names <- c("chromosome", "strand", "starts", "ends")
    msg <- paste0("e.g. exons_coordinates(gene1) <- list(chromosome =",
                " chr1, strand = '+', starts = c(10,50), end = c(30,80)).")
    object@exons_coordinates <- create_new_coordinates(value = value,
                                                    list_names = list_names,
                                                    msg = msg)
    validObject(object)
    return(object) })


#' @title Add/Remove Alternative Transcripts
#' @description Setter for alternative_transcripts attribute of `Gene` class
#' (and derived classes). By choosing an action, which can be 'add' or
#' 'remove', it's possible to add and remove, respectively, transcripts from
#' alternative_transcripts attribute.
#' @param object Object of class `Gene` (and derived classes).
#' @param action Action to be done: can be 'add' or 'remove'.
#' @param value Named list with transcript_id of the transcript to add/remove
#' and protein_coding, a `logical` to indicate the nature of the transcript.
#' @return The modified object.
#' @seealso
#' \link{alternative_transcripts<-,Gene-method}
#' @examples
#' # This is only the generic function, look in 'See Also' section
#' # for the method implementation.
#' cg1 <- CodingGene()
#' alternative_transcripts(cg1)
#' alternative_transcripts(cg1, action = 'add') <- list(transcript_id =
#' 'ENST12312312312', protein_coding = TRUE)
#' alternative_transcripts(cg1)
#' alternative_transcripts(cg1, action = 'remove') <- list(transcript_id =
#' 'ENST12312312312', protein_coding = TRUE)
#' alternative_transcripts(cg1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setGeneric("alternative_transcripts<-",
            function(object, action = 'add', value)
                standardGeneric("alternative_transcripts<-"))


utils::globalVariables("transcript_ensembl_id")

#' @title Add/Remove Alternative Transcripts
#' @description Setter for alternative_transcripts attribute of `Gene` class
#' (and derived classes). By choosing an action, which can be 'add' or
#' 'remove', it's possible to add and remove, respectively, transcripts from
#' alternative_transcripts attribute.
#' @param object Object of class `Gene` (and derived classes).
#' @param action Action to be done: can be 'add' or 'remove'.
#' @param value Named list with transcript_id of the transcript to add/remove
#' and protein_coding, a `logical` to indicate the nature of the transcript.
#' @return The modified object.
#' @examples
#' # Since Gene is a virtual class, the example is made using a derived
#' # class.
#' cg1 <- CodingGene()
#' alternative_transcripts(cg1)
#' alternative_transcripts(cg1, action = 'add') <- list(transcript_id =
#' 'ENST12312312312', protein_coding = TRUE)
#' alternative_transcripts(cg1)
#' alternative_transcripts(cg1, action = 'remove') <- list(transcript_id =
#' 'ENST12312312312', protein_coding = TRUE)
#' alternative_transcripts(cg1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("alternative_transcripts<-", "Gene",
            function(object, action, value){
    msg<-paste0("e.g. 'alternative_transcripts(gene1, action='add') <- list",
                    "(transcript_id = 'ENST12345678912',",
                    " protein_coding = TRUE)'.")
    if ( !action %in% c("add", "remove") || length(action) != 1) {
        err_msg <- paste0("Action argument accepts only 'add' and 'remove'",
                            " as input values, ", msg)
        stop(err_msg) }

    list_names <- c("transcript_id", "protein_coding")
    if (!all(names(value) %in% list_names)  || length(value) != 2) {
        err_msg <- paste0("The list element names, passed as arguments,",
                        " are < ",paste(list_names, collapse = ', '), " >, ",
                        msg)
        stop(err_msg) }

    row <- list(value$transcript_id, value$protein_coding)
    if (action == 'add') {
        df <- object@alternative_transcripts
        if (value$transcript_id %in% df$transcript_ensembl_id) {
            err_msg <- paste("Transcript already present among the",
                            "alternative transcripts of this gene.")
            stop(err_msg)}
        if (nrow(df)) {
            df <- rbind(df, row) }
        else {
            df <- data.frame(transcript_ensembl_id = row[[1]],
                            protein_coding = row[[2]]) }
    }
    else {
        df <- object@alternative_transcripts
        df <- dplyr::filter(df, transcript_ensembl_id != row[[1]])
    }

    df <- dplyr::filter(df, transcript_ensembl_id != 'ENST00000000000')
    object@alternative_transcripts <- unique(df)
    validObject(object)
    return(object) })
