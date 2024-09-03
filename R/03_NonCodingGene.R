#' @title NonCodingGene Virtual Class
#' @description This virtual class is the class from which
#' HousekeepingNonCodingGene and RegulatoryNonCodingGene classes are derived.
#' @slot pre_ncRNA_ensembl_id `character`. pre-ncRNA ENSEMBL ID.
#' @slot pre_ncRNA_coordinates `GRanges`.
#' pre-ncRNA genomic coordinates.
#' @return No objects can be created directly from this class since
#' it's virtual.
#' @seealso
#' \linkS4class{Gene}, \linkS4class{HousekeepingNonCodingGene},
#' \linkS4class{RegulatoryNonCodingGene}
#' @examples
#' # The NonCodingGene class is virtual, so no objects can be created
#' # directly from this class. It needs to be extended by other classes.
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setClass("NonCodingGene",
        contains = c("Gene", "VIRTUAL"),
        slots = c(
            pre_ncRNA_ensembl_id = "character",
            pre_ncRNA_coordinates = "GRanges"
        ),
        prototype = list(
            pre_ncRNA_ensembl_id = "ENST00000000000",
            pre_ncRNA_coordinates = GRanges()
        ))


#' @title Validity Check for NonCodingGene Virtual Class
#' @description The function checks that: each slot, except the special ones,
#' is of length 1 (or empty); pre-ncRNA ID is a valid ENSEMBL ID; all
#' coordinates are consistent w.r.t. each other and that gene sequence is
#' consistent with the provided coordinates. It also checks that no
#' protein-coding transcripts are present in alternative_transcripts slot.
#' @param object Object of class `NonCodingGene` (and derived classes).
#' @return `TRUE` if the object is valid, otherwise the function
#' throws an error via `stop()`.
#' @seealso
#' \link{Gene-validity}
#' @examples
#' # The NonCodingGene class is virtual, so no objects can be created
#' # directly from this class. It needs to be extended by other classes.
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @name NonCodingGene-validity

setValidity("NonCodingGene", function(object) {
    ### 1. Check on slots length (they must contain 1 element, except
    ###    the special ones).
    check_slots_length(obj = object, class_obj = "NonCodingGene",
                        mother_class_obj = "Gene")

    ### 2. Check on the correctness of precursor ncRNA ENSEMBL ID
    pre_ncRNA_id <- object@pre_ncRNA_ensembl_id
    check_ensembl_id(id = pre_ncRNA_id, prefix = 'ENST',
                    type = 'precursor ncRNA transcript')

    ### 3. Check that all coordinates share same chromosome and strand.
    check_coordinates_consistency(granges_list = list(
        gene_coordinates = object@gene_coordinates,
        exons_coordinates= object@exons_coordinates,
        pre_ncRNA_coordinates = object@pre_ncRNA_coordinates))

    ### 4. Check on precursor ncRNA coordinates w.r.t gene and exons
    ###    coordinates
    msg <- "pre-ncRNA coordinates can't exceed gene coordinates."
    check_boundaries(grange_in = object@pre_ncRNA_coordinates,
                    grange_out = object@gene_coordinates, msg = msg)

    check_in_exons(range_to_check = object@pre_ncRNA_coordinates,
                    exons = object@exons_coordinates,
                    type = "precursor non-coding RNA")

    ### 5. Check length gene sequence
    if (length(object@gene_sequence) &&
        !isEmpty(object@pre_ncRNA_coordinates)) {
        if (length(object@gene_sequence) <
            width(object@pre_ncRNA_coordinates)) {
            err_msg <- paste0("Gene sequence can't be shorter than",
                            " precursor non-coding RNA coordinates.")
            stop(err_msg) }
    }

    ### 6. Check on dataframe of slot 'alternative_transcripts'
    protein_coding_col <- object@alternative_transcripts$protein_coding
    if (any(protein_coding_col)) {
        err_msg <- paste0("A non-coding gene can't have a protein coding",
                        " transcript among its alternative transcripts.")
        stop(err_msg) }

    return(TRUE) })


# GETTERs for NonCodingGene

#' @title Get pre-ncRNA ENSEMBL ID
#' @description Getter for pre_ncRNA_ensembl_id attribute of `NonCodingGene`
#' class (and derived classes).
#' @param object Object of class `NonCodingGene` (and derived classes).
#' @return Value of pre_ncRNA_ensembl_id.
#' @seealso
#' \link{pre_ncRNA_id,NonCodingGene-method}
#' @examples
#' # This is only the generic function, look in 'See Also' section
#' # for the method implementation.
#' mir1 <- miRNAGene(pre_ncRNA_ensembl_id = 'ENST12312312300')
#' pre_ncRNA_id(mir1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setGeneric("pre_ncRNA_id",
            function(object) standardGeneric("pre_ncRNA_id"))


#' @title Get pre-ncRNA ENSEMBL ID
#' @description Getter for pre_ncRNA_ensembl_id attribute of `NonCodingGene`
#' class (and derived classes).
#' @param object Object of class `NonCodingGene` (and derived classes).
#' @return Value of pre_ncRNA_ensembl_id.
#' @examples
#' # Since NonCodingGene is a virtual class, the example is made using a
#' # derived class
#' mir1 <- miRNAGene(pre_ncRNA_ensembl_id = 'ENST12312312300')
#' pre_ncRNA_id(mir1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("pre_ncRNA_id", "NonCodingGene", function(object) {
    object@pre_ncRNA_ensembl_id })


#' @title Get pre-ncRNA Genomic Coordinates
#' @description Getter for pre_ncRNA_coordinates attribute of `NonCodingGene`
#' class (and derived classes).
#' @param object Object of class `NonCodingGene` (and derived classes).
#' @return `GRanges` object representing pre-ncRNA coordinates.
#' @seealso
#' \link{pre_ncRNA_coordinates,NonCodingGene-method}
#' @examples
#' # This is only the generic function, look in 'See Also' section
#' # for the method implementation.
#' mir1 <- miRNAGene(chromosome = 'chr1', strand = '+', pre_ncRNA_start=210,
#' pre_ncRNA_end = 800)
#' pre_ncRNA_coordinates(mir1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setGeneric("pre_ncRNA_coordinates",
            function(object) standardGeneric("pre_ncRNA_coordinates"))


#' @title Get pre-ncRNA Genomic Coordinates
#' @description Getter for pre_ncRNA_coordinates attribute of `NonCodingGene`
#' class (and derived classes).
#' @param object Object of class `NonCodingGene` (and derived classes).
#' @return `GRanges` object representing pre-ncRNA coordinates.
#' @examples
#' # Since NonCodingGene is a virtual class, the example is made using a
#' # derived class
#' mir1 <- miRNAGene(chromosome = 'chr1', strand = '+', pre_ncRNA_start=210,
#' pre_ncRNA_end = 800)
#' pre_ncRNA_coordinates(mir1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("pre_ncRNA_coordinates", "NonCodingGene", function(object) {
    object@pre_ncRNA_coordinates })


# SETTERs for NonCodingGene

#' @title Set pre-ncRNA ENSEMBL ID
#' @description Setter for pre_ncRNA_ensembl_id attribute of `NonCodingGene`
#' class (and derived classes).
#' @param object Object of class `NonCodingGene` (and derived classes).
#' @param value New value for pre_ncRNA_ensembl_id.
#' @return The modified object.
#' @seealso
#' \link{pre_ncRNA_id<-,NonCodingGene-method}
#' @examples
#' # This is only the generic function, look in 'See Also' section
#' # for the method implementation.
#' mir1 <- miRNAGene(pre_ncRNA_ensembl_id = 'ENST12312312300')
#' pre_ncRNA_id(mir1)
#' pre_ncRNA_id(mir1) <- 'ENST00099988812'
#' pre_ncRNA_id(mir1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setGeneric("pre_ncRNA_id<-",
            function(object, value) standardGeneric("pre_ncRNA_id<-"))


#' @title Set pre-ncRNA ENSEMBL ID
#' @description Setter for pre_ncRNA_ensembl_id attribute of `NonCodingGene`
#' class (and derived classes).
#' @param object Object of class `NonCodingGene` (and derived classes).
#' @param value New value for pre_ncRNA_ensembl_id.
#' @return The modified object.
#' @examples
#' # Since NonCodingGene is a virtual class, the example is made using a
#' # derived class
#' mir1 <- miRNAGene(pre_ncRNA_ensembl_id = 'ENST12312312300')
#' pre_ncRNA_id(mir1)
#' pre_ncRNA_id(mir1) <- 'ENST00099988812'
#' pre_ncRNA_id(mir1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("pre_ncRNA_id<-", "NonCodingGene", function(object, value) {
    object@pre_ncRNA_ensembl_id <- value
    validObject(object)
    return(object) })


#' @title Set pre-ncRNA Genomic Coordinates
#' @description Setter for pre_ncRNA_coordinates attribute of `NonCodingGene`
#' class (and derived classes).
#' @param object Object of class `NonCodingGene` (and derived classes).
#' @param value Named list with chromosome, strand, start and end values
#' to build the new pre-ncRNA coordinates.
#' @return The modified object.
#' @seealso
#' \link{pre_ncRNA_coordinates<-,NonCodingGene-method}
#' @examples
#' # This is only the generic function, look in 'See Also' section
#' # for the method implementation.
#' mir1 <- miRNAGene(chromosome = 'chr1', strand = '+', pre_ncRNA_start = 20,
#' pre_ncRNA_end = 300)
#' pre_ncRNA_coordinates(mir1)
#' pre_ncRNA_coordinates(mir1) <- list(chromosome = 'chr2', strand = '+',
#' start = 30, end = 50)
#' pre_ncRNA_coordinates(mir1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setGeneric("pre_ncRNA_coordinates<-",
            function(object, value)
                standardGeneric("pre_ncRNA_coordinates<-"))


#' @title Set pre-ncRNA Genomic Coordinates
#' @description Setter for pre_ncRNA_coordinates attribute of `NonCodingGene`
#' class (and derived classes).
#' @param object Object of class `NonCodingGene` (and derived classes).
#' @param value Named list with chromosome, strand, start and end values
#' to build the new pre-ncRNA coordinates.
#' @return The modified object.
#' @examples
#' # Since NonCodingGene is a virtual class, the example is made using a
#' # derived class
#' mir1 <- miRNAGene(chromosome = 'chr1', strand = '+', pre_ncRNA_start = 20,
#' pre_ncRNA_end = 300)
#' pre_ncRNA_coordinates(mir1)
#' pre_ncRNA_coordinates(mir1) <- list(chromosome = 'chr2', strand = '+',
#' start = 30, end = 50)
#' pre_ncRNA_coordinates(mir1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("pre_ncRNA_coordinates<-", "NonCodingGene", function(object, value){
    list_names <- c("chromosome", "strand", "start", "end")
    msg <- paste0(" e.g. 'pre_ncRNA_coordinates(gene1) <- list(chromosome =",
                    " chr1, strand = '+', start = 10, end = 50)'.")
    object@pre_ncRNA_coordinates <- create_new_coordinates(value = value,
                                                    list_names = list_names,
                                                    msg = msg)
    validObject(object)
    return(object) })
