#' @title rRNAGene Class
#' @description Class to represent genes encoding for ribosomal RNAs.
#' @slot rRNA_mature_sequence `RNAString`. rRNA mature sequence.
#' @slot rRNA_type `character`. Ribosomal RNA type (e.g. 5S) .
#' @return An object of class `rRNAGene`.
#' @seealso
#' \linkS4class{HousekeepingNonCodingGene}
#' @examples
#' # An empty rRNAGene object
#' r1 <- new("rRNAGene")
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setClass("rRNAGene",
        contains = "HousekeepingNonCodingGene",
        slots = c(
            rRNA_mature_sequence = "RNAString",
            rRNA_type = "character"
        ),
        prototype = list(
            rRNA_mature_sequence = RNAString(),
            rRNA_type = NA_character_
        ) )


#' @title Validity Check for rRNAGene Class
#' @description The function checks that: each slot, except the special ones,
#' is of length 1 (or empty); maure tRNA sequence is compatible with both the
#' defined coordinates and the gene sequence and finally, that rRNA type and
#' chromosome are consistent (e.g. rRNA 5S genes are only on chromosome 1).
#' @param object Object of class `rRNAGene`.
#' @return `TRUE` if the object is valid, otherwise the function
#' throws an error via `stop()`.
#' @seealso
#' \link{HousekeepingNonCodingGene-validity}
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @name rRNAGene-validity

setValidity("rRNAGene", function(object) {
    ### 1. Check on slots length (they must contain 1 element, except
    ###    the special ones).
    check_slots_length(obj = object, class_obj = "rRNAGene",
                        mother_class_obj = "HousekeepingNonCodingGene",
                        special_slots = "rRNA_mature_sequence")
    ### 2. Check length of mature rRNA sequence w.r.t
    ###    gene/exons/ncRNA precursor
    check_mature_length(object=object,mature_sequence='rRNA_mature_sequence')
    ### 3. Check on rRNA mature sequence length if coordinates are
    ### not defined but gene sequence is.
    if (length(object@gene_sequence)) {
        if (length(object@gene_sequence) <
            length(object@rRNA_mature_sequence)) {
            err_msg <- paste0("The sequence of the gene can't be",
                            " shorter than its mature rRNA transcript.")
            stop(err_msg) }
    }
    ### 4. Check validity of rRNA_type.
    human_rRNA_types <- c("5S","5.8S","18S","28S","12S","16S")
    if (!object@rRNA_type %in% c(human_rRNA_types, NA)) {
        err_msg <- paste0("Please provide a valid human rRNA type.",
                    " Choose among {", paste(human_rRNA_types,
                                            collapse = ', '), "}")
        stop(err_msg) }
    ### 5. Check consistency of rRNA type with chromosome.
    coordinates <- list(object@gene_coordinates,
                        object@exons_coordinates,
                        object@pre_ncRNA_coordinates)
    chromosome <- check_coordinates_consistency(coordinates,
                                                return_chr = TRUE)

    if (length(chromosome) && !is.na(object@rRNA_type)) {
        if ( object@rRNA_type %in% human_rRNA_types[c(5,6)]
            && (chromosome != "chrM") ) {
                err_msg <- paste0("This rRNA type is encoded by",
                            " mitochondrial chromosome,",
                            " not nuclear chromosomes.")
                stop(err_msg) }
        if (object@rRNA_type %in% human_rRNA_types[c(2,3,4)]
            && (!chromosome %in% c("chr13","chr14","chr15","chr21","chr22"))){
                err_msg <- paste0("This rRNA type is encoded only by genes",
                        " on chromosomes 13, 14, 15, 21 and 22.")
                stop(err_msg) }
        if (object@rRNA_type == "5S" && chromosome != "chr1") {
            err_msg <- paste0("This rRNA type is encoded only by genes on",
                            " chromosome 1.")
            stop(err_msg)} }
    return(TRUE) })


# GETTERs for rRNA Gene

#' @title Get Mature rRNA Sequence
#' @description Getter for rRNA_mature_sequence attribute of `rRNAGene` class.
#' @param object Object of class `rRNAGene`.
#' @return Mature rRNA sequence.
#' @examples
#' r1 <- rRNAGene(rRNA_mature_sequence = 'AUGCA')
#' mature_sequence(r1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("mature_sequence", "rRNAGene", function(object) {
    object@rRNA_mature_sequence })

#' @title Get rRNA Type
#' @description Getter for rRNA_type attribute of `rRNAGene` class.
#' @param object Object of class `rRNAGene`.
#' @return rRNA type value.
#' @seealso
#' \link{rRNA_type,rRNAGene-method}
#' @examples
#' # This is only the generic function, look in 'See Also' section
#' # for the method implementation.
#' r1 <- rRNAGene(rRNA_type='5S')
#' rRNA_type(r1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setGeneric("rRNA_type", function(object) standardGeneric("rRNA_type"))


#' @title Get rRNA Type
#' @description Getter for rRNA_type attribute of `rRNAGene` class.
#' @param object Object of class `rRNAGene`.
#' @return rRNA type value.
#' @examples
#' r1 <- rRNAGene(rRNA_type = '5S')
#' rRNA_type(r1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("rRNA_type", "rRNAGene", function(object) {
    object@rRNA_type })


# SETTERs for rRNA Gene

#' @title Set Mature rRNA Sequence
#' @description Setter for rRNA_mature_sequence attribute of `rRNAGene` class.
#' @param object Object of class `rRNAGene`.
#' @param value New value for rRNA_mature_sequence.
#' @return The modified object.
#' @examples
#' r1 <- rRNAGene(rRNA_mature_sequence = 'AUUA')
#' mature_sequence(r1)
#' mature_sequence(r1) <- 'AAAGCA'
#' mature_sequence(r1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("mature_sequence<-", "rRNAGene", function(object, value) {
    object@rRNA_mature_sequence <- RNAString(as.character(value))
    validObject(object)
    return(object) })


#' @title Set rRNA Type
#' @description Setter for rRNA_type attribute of `rRNAGene` class.
#' @param object Object of class `rRNAGene`.
#' @param value New value for rRNA_type
#' @return The modified object.
#' @seealso
#' \link{rRNA_type<-,rRNAGene-method}
#' @examples
#' # This is only the generic function, look in 'See Also' section
#' # for the method implementation.
#' r1 <- rRNAGene(rRNA_type = '5S')
#' rRNA_type(r1)
#' rRNA_type(r1) <- '18S'
#' rRNA_type(r1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setGeneric("rRNA_type<-",
            function(object, value) standardGeneric("rRNA_type<-"))


#' @title Set rRNA Type
#' @description Setter for rRNA_type attribute of `rRNAGene` class.
#' @param object Object of class `rRNAGene`.
#' @param value New value for rRNA_type
#' @return The modified object.
#' @examples
#' r1 <- rRNAGene(rRNA_type = '5S')
#' rRNA_type(r1)
#' rRNA_type(r1) <- '18S'
#' rRNA_type(r1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("rRNA_type<-", "rRNAGene", function(object, value) {
    object@rRNA_type <- value
    validObject(object)
    return(object) })


#' @title Create an Object of Class rRNAGene
#' @usage
#' rRNAGene(
#'     gene_ensembl_id = "ENSG00000000000",
#'     hugo_symbol = NA_character_,
#'     gene_complete_name = NA_character_,
#'     gene_description = NA_character_,
#'     strand = NA,
#'     chromosome = NA,
#'     gene_start = NA,
#'     gene_end = NA,
#'     gene_sequence = DNAString(),
#'     exons_starts = NA,
#'     exons_ends = NA,
#'     pre_ncRNA_ensembl_id = "ENST00000000000",
#'     pre_ncRNA_start = NA,
#'     pre_ncRNA_end = NA,
#'     essentiality_score = NA_real_,
#'     ubiquitous_expression = NA,
#'     rRNA_mature_sequence = RNAString(),
#'     rRNA_type = NA_character_ )
#' @description Constructor function for `rRNAGene` class.
#' @param gene_ensembl_id Gene ENSEMBL ID.
#' @param hugo_symbol Gene Hugo symbol.
#' @param gene_complete_name Gene complete name.
#' @param gene_description Gene description.
#' @param strand Strand on which the gene is annotated.
#' @param chromosome Chromosome on which the gene is annotated.
#' @param gene_start Gene start coordinate.
#' @param gene_end Gene end coordinate.
#' @param gene_sequence Gene sequence.
#' @param exons_starts Exons starts coordinates.
#' @param exons_ends Exons ends coordinates.
#' @param pre_ncRNA_ensembl_id pre-ncRNA ENSEMBL ID.
#' @param pre_ncRNA_start pre-ncRNA start coordinate.
#' @param pre_ncRNA_end pre-ncRNA end coordinate.
#' @param essentiality_score Gene essentiality score.
#' @param ubiquitous_expression `logical` for ubiquitous expression.
#' @param rRNA_mature_sequence rRNA mature sequence.
#' @param rRNA_type Type of rRNA.
#' @return New rRNAGene object.
#' @examples
#' # An empty rRNAGene object
#' r1 <- rRNAGene()
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

rRNAGene <- function(gene_ensembl_id = "ENSG00000000000",
                    hugo_symbol = NA_character_,
                    gene_complete_name = NA_character_,
                    gene_description = NA_character_,
                    strand = NA, chromosome = NA,
                    gene_start = NA, gene_end = NA,
                    gene_sequence = DNAString(),
                    exons_starts = NA, exons_ends = NA,
                    pre_ncRNA_ensembl_id = "ENST00000000000",
                    pre_ncRNA_start = NA, pre_ncRNA_end = NA,
                    essentiality_score = NA_real_,
                    ubiquitous_expression = NA,
                    rRNA_mature_sequence = RNAString(),
                    rRNA_type = NA_character_) {

    coordinates <- constructor_helper(chrom = chromosome, str = strand,
                                    g_st = gene_start, g_en = gene_end,
                                    ex_st = exons_starts, ex_en=exons_ends,
                                    pre_st = pre_ncRNA_start,
                                    pre_en = pre_ncRNA_end,
                                    coding = FALSE)
    new("rRNAGene",
        gene_ensembl_id = gene_ensembl_id,
        hugo_symbol = hugo_symbol,
        gene_complete_name = gene_complete_name,
        gene_description = gene_description,
        gene_coordinates = coordinates[[1]],
        gene_sequence = DNAString(gene_sequence),
        exons_coordinates = coordinates[[2]],
        pre_ncRNA_ensembl_id = pre_ncRNA_ensembl_id,
        pre_ncRNA_coordinates = coordinates[[3]],
        essentiality_score = essentiality_score,
        ubiquitous_expression = ubiquitous_expression,
        rRNA_mature_sequence = RNAString(rRNA_mature_sequence),
        rRNA_type = rRNA_type) }

